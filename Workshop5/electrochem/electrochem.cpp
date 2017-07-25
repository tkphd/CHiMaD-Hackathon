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
	Field 1 is electrostatic potential.
*/

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
		vector<double> gradc = gradient(GRID, x, 0);

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

template<int dim,typename T>
unsigned int RedBlackGaussSeidel(const grid<dim,vector<T> >& oldGrid, const double timestep, grid<dim,vector<T> >& newGrid)
{
	double gridSize = static_cast<double>(nodes(oldGrid));
	#ifdef MPI_VERSION
	double localGridSize = gridSize;
	MPI::COMM_WORLD.Allreduce(&localGridSize, &gridSize, 1, MPI_DOUBLE, MPI_SUM);
	#endif

	double dV = 1.0;
	double lapWeight = 0.0;
	for (int d=0; d<dim; d++) {
		dV *= dx(oldGrid,d);
		lapWeight += 2.0 / std::pow(dx(oldGrid,d),2.0); // dim=2 -> lapWeight = 4/h^2 if dy=dx=h
	}

	double residual = 1.0;
	unsigned int iter = 0;

	while (iter<max_iter && residual>tolerance) {
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
			#pragma omp parallel for schedule(dynamic)
			for (int n=0; n<nodes(oldGrid); n++) {
				vector<int> x = position(oldGrid,n);

				// Determine tile color
				int x_sum=0;
				for (int d=0; d<dim; d++)
					x_sum += x[d];
				if (x_sum%2 == color)
					continue;

				const T& cOld = oldGrid(n)[cid];
				const T& cGuess = newGrid(n)[cid]; // value from last "guess" iteration
				const T& uGuess = newGrid(n)[uid];
				T pGuess = newGrid(n)[pid];

				if (x[0] == g1(oldGrid, 0) - 1)
					pGuess = std::sin(dx(oldGrid, 1)/7.0  * x[1]);
				else if (x[0] == g0(oldGrid, 0))
					pGuess = 0.0;

				// A is defined by the last guess, stored in newGrid(n). It is a 3x3 matrix.
				const double A12 = lapWeight * timestep * M;

				const double A21 = -kappa * lapWeight - dfcontractivedc(cGuess, 1.0);

				const double A31 = k / epsilon;
				const double A33 = -lapWeight;

				// B is defined by the last value, stored in oldGrid(n), and the last guess, stored in newGrid(n). It is a 3x1 column.
				const double B1 = cOld + timestep * M * fringe_laplacian(newGrid, x, uid);
				const double B2 = dfexpansivedc(cGuess, pGuess) - kappa * fringe_laplacian(newGrid, x, cid);
				const double B3 = -fringe_laplacian(newGrid, x, pid);

				// Solve the iteration system AX=B using Cramer's rule
				const double detA = A33 - A12*A21*A33;

				const double detA1 = B1*A33 - A12*B2*A33;
				const T cNew = (std::fabs(detA) > tolerance) ? detA1 / detA : 0.0;

				const double detA2 = B2*A33 - B1*A21*A33;
				const T uNew = (std::fabs(detA) > tolerance) ? detA2 / detA : 0.0;

				T pNew = 0.0;

				if (x[0] == g1(oldGrid, 0) - 1)
					pNew = std::sin(dx(oldGrid, 1)/7.0  * x[1]);
				else if (x[0] == g0(oldGrid, 0))
					pNew = 0.0;
				else {
					const double detA3 = B3 + A12*B2*A31 - B1*A31 - A12*A21*B3;
					pNew = (std::fabs(detA) > tolerance) ? detA3 / detA : 0.0;
				}

				// (Don't) Apply relaxation
				newGrid(n)[cid] = omega*cNew + (1.0 - omega)*cGuess;
				newGrid(n)[uid] = omega*uNew + (1.0 - omega)*uGuess;
				newGrid(n)[pid] = omega*pNew + (1.0 - omega)*pGuess;

			}
			ghostswap(newGrid);   // fill in the ghost cells; does nothing in serial
		}

		iter++;

		/*  ==== RESIDUAL ====
		    The residual is computed from the original matrix form, Ax=b:
		    any Old term goes into B, while any New term goes in AX. Note that
		    this is not the iteration matrix, it is the original system of equations.
		*/

		if (iter<residual_step || iter%residual_step==0) {
			double normB = 0.0;
			residual = 0.0;

			#pragma omp parallel for schedule(dynamic)
			for (int n=0; n<nodes(oldGrid); n++) {
				vector<int> x = position(oldGrid,n);
				vector<T> lap = laplacian(newGrid, x);

				const T& cOld = oldGrid(n)[cid];
				const T& cNew = newGrid(n)[cid];
				const T& uNew = newGrid(n)[uid];
				const T& pNew = newGrid(n)[pid];

				// Plug iteration results into original system of equations
				const double AX1 = cNew;
				const double AX2 = uNew;
				const double AX3 = lap[pid];

				const double B1 = cOld + timestep * M * lap[uid];
				const double B2 = dfchemdc(cNew) + 2.0*dfelecdc(pNew) - kappa*lap[cid];
				const double B3 = -k / epsilon * cNew;

				// Compute the Error from parts of the solution
				const double R1 = B1 - AX1;
				const double R2 = B2 - AX2;
				const double R3 = B3 - AX3;

				const double error = R1*R1 + R2*R2 + R3*R3;

				#pragma omp critical
				{
					residual += error;
					normB += B1*B1 + B2*B2 + B3*B3;
				}
			}

			#ifdef MPI_VERSION
			double localResidual=residual;
			MPI::COMM_WORLD.Allreduce(&localResidual, &residual, 1, MPI_DOUBLE, MPI_SUM);
			double localNormB=normB;
			MPI::COMM_WORLD.Allreduce(&localNormB, &normB, 1, MPI_DOUBLE, MPI_SUM);
			#endif

			residual = (std::fabs(normB) > tolerance)? sqrt(residual/normB)/(3.0*gridSize) : 0.0;
		}
	}

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
		GRID2D initGrid(3,0,100,0,100);
		for (int d=0; d<dim; d++) {
			// Set grid resolution
			dx(initGrid,d) = deltaX;

			// Set Neumann (zero-flux) boundary conditions
			if (x0(initGrid, d) == g0(initGrid, d))
				b0(initGrid, d) = Neumann;
			if (x1(initGrid, d) == g1(initGrid, d))
				b1(initGrid, d) = Neumann;
		}

		for (int n=0; n<nodes(initGrid); n++) {
			vector<int> x = position(initGrid, n);
			// composition field
			initGrid(n)[cid] = cheminit(dx(initGrid,0) * x[0], dx(initGrid,1) * x[1]);

			// charge field
			if (x[0] == g1(initGrid, 0) - 1)
				initGrid(n)[pid] = std::sin(dx(initGrid, 1)/7.0  * x[1]);
			else
				initGrid(n)[pid] = 0.0;
		}

		ghostswap(initGrid);

		for (int n=0; n<nodes(initGrid); n++) {
			// chemical potential field
			vector<int> x = position(initGrid, n);
			const double& c = initGrid(n)[cid];
			const double& p = initGrid(n)[pid];
			initGrid(n)[uid] = dfchemdc(c) + 2.0 * dfelecdc(p) - kappa * laplacian(initGrid, x, cid);
		}

		output(initGrid,filename);

		double F = Helmholtz(initGrid);

		if (rank == 0)
			std::cout << F << std::endl;
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

	grid<dim,vector<T> > newGrid(oldGrid);

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

	for (int step=0; step<steps; step++) {
		if (rank==0)
			print_progress(step, steps);

		ghostswap(oldGrid);
		ghostswap(newGrid);

		unsigned int iter = RedBlackGaussSeidel(oldGrid, dt, newGrid);

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
	}

	ghostswap(oldGrid);

	double F = Helmholtz(oldGrid);

	if (rank == 0)
		std::cout << F << std::endl;
}

} // MMSP
#endif

#include"MMSP.main.hpp"
