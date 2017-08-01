// electrochem.cpp
// Algorithms for 2D and 3D electrochem model
// Questions/comments to trevor.keller@nist.gov (Trevor Keller)

#ifndef CAHNHILLIARD_UPDATE
#define CAHNHILLIARD_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"electrochem.hpp"
#include"energy.hpp"

namespace MMSP {

/*
	Field 0 is composition.
	Field 1 is chemical potential.
	Field 2 is electrostatic potential.
*/

// Numerical parameters
const double deltaX = 0.2;
const double dt = 0.05;
const double CFL = (24.0*M*kappa) / std::pow(deltaX, 4);

template <int dim,typename T>
double Helmholtz(const grid<dim,vector<T> >& GRID)
{
	double dV = 1.0;
	for (int d=0; d<dim; d++)
		dV *= dx(GRID, d);

	double fchem = 0.0;
	double felec = 0.0;
	double fgrad = 0.0;

	for (int n=0; n<nodes(GRID); n++) {
		vector<int> x = position(GRID, n);
		vector<double> gradc = gradient(GRID, x, cid);

		fchem += chemenergy(GRID(n)[cid]);
		felec += elecenergy(GRID(n)[cid], GRID(n)[pid]);
		fgrad += gradc * gradc;
	}

	double F = dV*(fchem + felec + 0.5*kappa*fgrad); // Equation 5

	#ifdef MPI_VERSION
	double myF(F);
	MPI::COMM_WORLD.Allreduce(&myF, &F, 1, MPI_DOUBLE, MPI_SUM);
	#endif

	return F;
}

// Discrete Laplacian operator missing the central value, for implicit source terms
template<int dim, typename T>
double fringe_laplacian(const MMSP::grid<dim,MMSP::vector<T> >& GRID, const MMSP::vector<int>& x, const int field)
{
	double laplacian = 0.0;
	MMSP::vector<int> s = x;

	for (int i=0; i<dim; i++) {
		double weight = 1.0 / pow(dx(GRID, i), 2.0);

		s[i] += 1;
		const T& yh = GRID(s)[field];
		s[i] -= 2;
		const T& yl = GRID(s)[field];
		s[i] += 1;

		laplacian += weight * (yh + yl);
	}

	return laplacian;
}

template<int dim,typename T>
bool pointIsOnRightBoundary(const MMSP::grid<dim,T> GRID, const MMSP::vector<int>& x)
{
	return (x[0] == g1(GRID, 0) - 1);
}

template<int dim,typename T>
bool pointIsOnLeftBoundary(const MMSP::grid<dim,T>& GRID, const MMSP::vector<int>& x)
{
	return (x[0] == g0(GRID, 0));
}

template<int dim,typename T>
unsigned int RedBlackGaussSeidel(const grid<dim,vector<T> >& oldGrid, grid<dim,vector<T> >& newGrid)
{
	#ifdef DEBUG
	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

	std::ofstream of;
	if (rank == 0)
		of.open("iter.log", std::ofstream::out | std::ofstream::app); // new results will be appended
	#endif

	double gridSize = static_cast<double>(nodes(oldGrid));
	#ifdef MPI_VERSION
	double localGridSize = gridSize;
	MPI::COMM_WORLD.Allreduce(&localGridSize, &gridSize, 1, MPI_DOUBLE, MPI_SUM);
	#endif

	double dV = 1.0;
	double lapWeight = 0.0;
	for (int d=0; d<dim; d++) {
		dV *= dx(oldGrid,d);
		lapWeight += 2.0 / std::pow(dx(oldGrid, d), 2.0);
	}

	double residual = 2.0;
	unsigned int iter = 0;

	while (iter < max_iter && residual > tolerance) {
		/*  ==== RED-BLACK GAUSS SEIDEL ====
		    Iterate over a checkerboard, updating first red then black tiles.
		    This method eliminates the third "guess" grid, and should converge faster.
		    In 2-D and 3-D, if the sum of indices is even, then the tile is Red; else, Black.

			This method solves the linear system of equations,
		    /  1  a12  0  \ / x1 \   / b1 \
		    | a21  1   0  | | x2 | = | b2 |
			\ a31  0  a33 / \ x3 /   \ b3 /
		*/

		for (int color=1; color>-1; color--) {
			// If color==1, skip BLACK tiles, which have Σx[d] odd
			// If color==0, skip RED tiles, which have Σx[d] even
			#ifdef _OPENMP
			#pragma omp parallel for
			#endif
			for (int n=0; n<nodes(oldGrid); n++) {
				vector<int> x = position(oldGrid,n);

				// Determine tile color
				int x_sum=0;
				for (int d=0; d<dim; d++)
					x_sum += x[d];
				if (x_sum % 2 == color)
					continue;

				const T cOld = oldGrid(n)[cid];
				const T pOld = oldGrid(n)[pid];

				const T cGuess = newGrid(n)[cid]; // value from last "guess" iteration
				const T uGuess = newGrid(n)[uid];
				      T pGuess = newGrid(n)[pid];

				if (pointIsOnRightBoundary(oldGrid, x))
					pGuess = std::sin(dx(oldGrid, 1)/7.0  * x[1]);
				else if (pointIsOnLeftBoundary(oldGrid, x))
					pGuess = 0.0;

				// A is defined by the last guess, stored in newGrid(n). It is a 3x3 matrix.
				const double a12 = lapWeight * dt * M;
				const double a21 = -kappa * lapWeight - dfcontractivedc(cGuess, 1.0);
				const double a33 = -lapWeight;

				// B is defined by the last value, stored in oldGrid(n), and the last guess, stored in newGrid(n). It is a 3x1 column.
				const double flapC = fringe_laplacian(newGrid, x, cid);
				const double flapU = fringe_laplacian(newGrid, x, uid);
				const double flapP = fringe_laplacian(newGrid, x, pid);

				const double b1 = cOld + dt * M * flapU;
				const double b2 = dfexpansivedc(cOld, pOld) - kappa * flapC;
				const double b3 = -k * cOld / epsilon - flapP;

				// Solve the iteration system AX=B using Cramer's rule
				const double detA  = a33 * (1. - a12 * a21);
				const double detA1 = a33 * (b1 - a12 * b2 );
				const double detA2 = a33 * (b2 - b1  * a21);
				const double detA3 = b3  * (1. - a12 * a21);
				                // - a31 * (b1 - a12 * b2 );

				const T cNew = detA1 / detA;
				const T uNew = detA2 / detA;
				      T pNew = detA3 / detA;

				if (pointIsOnRightBoundary(newGrid, x))
					pNew = pGuess;
				else if (pointIsOnLeftBoundary(newGrid, x))
					pNew = 0.0;

				// (Don't) Apply relaxation
				newGrid(n)[cid] = omega * cNew + (1.0 - omega) * cGuess;
				newGrid(n)[uid] = omega * uNew + (1.0 - omega) * uGuess;
				newGrid(n)[pid] = omega * pNew + (1.0 - omega) * pGuess;

			}
			ghostswap(newGrid);   // fill in the ghost cells; does nothing in serial
		}

		iter++;

		/*  ==== RESIDUAL ====
		    The residual is computed from the original matrix form, Ax=b:
		    any Old term goes into B, while any New term goes in AX. Note that
		    this is not the iteration matrix, it is the original system of equations.
		*/

		if (iter < residual_step || iter % residual_step == 0) {
			double normB = 0.0;
			residual = 0.0;

			#ifdef _OPENMP
			#pragma omp parallel for
			#endif
			for (int n=0; n<nodes(newGrid); n++) {
				vector<int> x = position(newGrid,n);
				vector<T> lap = laplacian(newGrid, x);

				const T cOld = oldGrid(n)[cid];
				const T pOld = oldGrid(n)[pid];

				const T cNew = newGrid(n)[cid];
				const T uNew = newGrid(n)[uid];

				// Plug iteration results into original system of equations
				const double Ax1 = cNew - dt * M * lap[uid];
				const double Ax2 = uNew - dfcontractivedc(cNew, cNew) + kappa * lap[cid];
				const double Ax3 = lap[pid];

				const double b1 = cOld;
				const double b2 = dfexpansivedc(cOld, pOld);
				const double b3 = -k * cOld / epsilon;

				// Compute the Error from parts of the solution
				const double r1 = b1 - Ax1;
				const double r2 = b2 - Ax2;
				const double r3 = (pointIsOnLeftBoundary(newGrid, x) || pointIsOnRightBoundary(newGrid, x)) ? 0.0 : b3 - Ax3;

				const double error  = r1*r1 + r2*r2 + r3*r3;
				const double source = b1*b1 + b2*b2 + b3*b3;

				#ifdef _OPENMP
				#pragma omp atomic
				#endif
				residual += error;

				#ifdef _OPENMP
				#pragma omp atomic
				#endif
				normB += source;
			}

			#ifdef MPI_VERSION
			double localResidual = residual;
			MPI::COMM_WORLD.Allreduce(&localResidual, &residual, 1, MPI_DOUBLE, MPI_SUM);
			double localNormB = normB;
			MPI::COMM_WORLD.Allreduce(&localNormB, &normB, 1, MPI_DOUBLE, MPI_SUM);
			#endif

			residual = (std::fabs(normB) > tolerance) ? sqrt(residual/normB)/(3.0*gridSize) : 0.0;

			#ifdef DEBUG
			double F = Helmholtz(oldGrid);

			if (rank == 0)
				of << iter << '\t' << residual << '\t' << F << '\n';
			#endif
		}

	}

	#ifdef DEBUG
	if (rank == 0)
		of.close();
	#endif

	return iter;
}

void generate(int dim, const char* filename)
{
	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

	if (dim != 2 && rank == 0) {
		std::cerr << "ERROR: CHiMaD problems are 2-D, only!" << std::endl;
		std::exit(-1);
	}

	/*	Grid contains three fields:
		0: composition
		1: chemical potential
		2: electrostatic potential  */

	if (dim==2) {
		const int L = 100 / deltaX;
		GRID2D initGrid(3, 0,L, 0,L);
		for (int d=0; d<dim; d++) {
			// Set grid resolution
			dx(initGrid,d) = deltaX;

			// Set Neumann (zero-flux) boundary conditions
			if (x0(initGrid, d) == g0(initGrid, d))
				b0(initGrid, d) = Neumann;
			if (x1(initGrid, d) == g1(initGrid, d))
				b1(initGrid, d) = Neumann;
		}

		if (rank == 0)
			std::cout << "Timestep is " << dt << ". CFL is " << CFL << ". Run " << 1.0 / dt << " per unit time." << std::endl;

		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for (int n=0; n<nodes(initGrid); n++) {
			vector<int> x = position(initGrid, n);
			// composition field
			initGrid(n)[cid] = cheminit(dx(initGrid,0) * x[0], dx(initGrid,1) * x[1]);

			// charge field
			if (pointIsOnRightBoundary(initGrid, x))
				initGrid(n)[pid] = std::sin(dx(initGrid, 1)/7.0  * x[1]);
			else
				initGrid(n)[pid] = 0.0;
		}

		ghostswap(initGrid);

		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for (int n=0; n<nodes(initGrid); n++) {
			// chemical potential field
			vector<int> x = position(initGrid, n);
			const double& c = initGrid(n)[cid];
			const double& p = initGrid(n)[pid];
			const double lapC = laplacian(initGrid, x, cid);
			initGrid(n)[uid] = dfcontractivedc(c, c) + dfexpansivedc(c, p) - kappa * lapC;
		}

		output(initGrid,filename);

		double F = Helmholtz(initGrid);
		std::ofstream of;

		#ifdef DEBUG
		if (rank == 0) {
			of.open("iter.log");
			of << "iter\tres\tF\n";
			of << 0 << '\t' << 0 << '\t' << F << '\n';
			of.close();
		}
		#else
		if (rank == 0) {
			of.open("energy_dx02.tsv");
			of << "t\tF\n";
			of << 0 << '\t' << F << '\n';
			of.close();
		}
		#endif
	}
}

template <int dim, typename T>
void update(grid<dim,vector<T> >& oldGrid, int steps)
{
	/*	Grid contains three fields:
		0: composition
		1: chemical potential
		2: electrostatic potential  */

	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

	#ifndef DEBUG
	static double elapsed = 0.0;
	#endif

	ghostswap(oldGrid);

	grid<dim,vector<T> > newGrid(oldGrid);
	newGrid.copy(oldGrid);

	// Make sure the grid spacing is correct
	for (int d=0; d<dim; d++) {
		dx(oldGrid,d) = deltaX;
		dx(newGrid,d) = deltaX;

		// Set Neumann (zero-flux) boundary conditions
		if (x0(oldGrid, d) == g0(oldGrid, d)) {
			b0(oldGrid, d) = Neumann;
			b0(newGrid, d) = Neumann;
		}
		if (x1(oldGrid, d) == g1(oldGrid, d)) {
			b1(oldGrid, d) = Neumann;
			b1(newGrid, d) = Neumann;
		}
	}

	#ifndef DEBUG
	std::ofstream of;
	if (rank == 0)
		of.open("energy_dx02.tsv", std::ofstream::out | std::ofstream::app); // new results will be appended
	#endif

	for (int step=0; step<steps; step++) {
		if (rank==0)
			print_progress(step, steps);

		ghostswap(newGrid);

		unsigned int iter = RedBlackGaussSeidel(oldGrid, newGrid);

		#ifdef MPI_VERSION
		unsigned int myit(iter);
		MPI::COMM_WORLD.Allreduce(&myit, &iter, 1, MPI_UNSIGNED, MPI_MAX);
		#endif

		if (iter >= max_iter) {
			if (rank==0)
				std::cerr << "Solver stagnated on step " << step << ". Aborting." << std::endl;
				MMSP::Abort(-1);
		}

		swap(oldGrid, newGrid);
		ghostswap(oldGrid);

		#ifndef DEBUG
		elapsed += dt;
		const double F = Helmholtz(oldGrid);
		if (rank == 0)
			of << elapsed << '\t' << F << '\n';
		#endif
	}

	#ifndef DEBUG
	if (rank == 0)
		of.close();
	#endif

}

} // MMSP

#endif

#include"MMSP.main.hpp"
