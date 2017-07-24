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

		fchem += chemenergy(GRID(n)[0]);
		felec += elecenergy(GRID(n)[0], GRID(n)[1]);
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
unsigned int RedBlackGaussSeidel(const grid<dim,vector<T> >& oldGrid, grid<dim,vector<T> >& newGrid)
{
	double gridSize = static_cast<double>(nodes(oldGrid));
	#ifdef MPI_VERSION
	double localGridSize = gridSize;
	MPI::COMM_WORLD.Allreduce(&localGridSize, &gridSize, 1, MPI_DOUBLE, MPI_SUM);
	#endif

	double dV = 1.0;
	for (int d=0; d<dim; d++)
		dV *= dx(oldGrid,d);

	double residual = 1.0;
	unsigned int iter = 0;

	while (iter<max_iter && residual>tolerance) {

		/*  ==== RED-BLACK GAUSS SEIDEL ====
		    Iterate over a checkerboard, updating first red then black tiles.
		    This method eliminates the third "guess" grid, and should converge faster.
		    In 2-D and 3-D, if the sum of indices is even, then the tile is Red; else, Black.
		*/

		for (int color = 1; color > -1; color--) {
			// If color==1, skip BLACK tiles, which have Σx[d] odd
			// If color==0, skip RED tiles, which have Σx[d] even
			#pragma omp parallel for schedule(dynamic)
			for (int n=0; n<nodes(oldGrid); n++) {
				vector<int> x = position(oldGrid,n);
				int x_sum = 0;
				for (int d = 0; d < dim; d++)
					x_sum += x[d];
				if (x_sum % 2 == color)
					continue;

				if (x[0] == g1(oldGrid, 0) - 1)
					newGrid(n)[1] = std::sin(dx(oldGrid, 1)/7.0  * x[1]);
				else if (x[0] == g0(oldGrid, 0))
					newGrid(n)[1] = 0.0;
				else {
					const T& cOld = oldGrid(n)[0];
					const T pGuess = newGrid(n)[1]; // value from last "guess" iteration

					const T pNew = (dV / 2.0*dim) * (fringe_laplacian(newGrid, x, 1) + k * cOld / epsilon);

					// Apply over-relaxation
					newGrid(n)[1] = omega * pNew + (1.0 - omega) * pGuess;
				}
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
				const T lap = laplacian(newGrid, x, 1);
				const T cOld = oldGrid(n)[0];

				const T AX = lap;
				const T B = -k * cOld / epsilon;
				const double R = B - AX;
				const double error = R * R;

				#pragma omp critical
				{
					residual += error;
					normB += B * B;
				}
			}

			#ifdef MPI_VERSION
			double localResidual=residual;
			MPI::COMM_WORLD.Allreduce(&localResidual, &residual, 1, MPI_DOUBLE, MPI_SUM);
			double localNormB=normB;
			MPI::COMM_WORLD.Allreduce(&localNormB, &normB, 1, MPI_DOUBLE, MPI_SUM);
			#endif

			residual = sqrt(residual/normB) / gridSize;
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
		std::cerr<<"ERROR: CHiMaD problems are 2-D, only!"<<std::endl;
		std::exit(-1);
	}

	if (dim==2) {
		GRID2D initGrid(2,0,100,0,100);
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
			initGrid(n)[0] = cheminit(dx(initGrid,0)*x[0], dx(initGrid,1)*x[1]);
			if (x[0] == g1(initGrid, 0) - 1)
				initGrid(n)[1] = std::sin(dx(initGrid, 1)/7.0  * x[1]);
			else
				initGrid(n)[1] = 0.0;
		}

		unsigned int iter = RedBlackGaussSeidel(initGrid, initGrid);
		if (iter >= max_iter) {
			if (rank==0)
				std::cerr << "Solver stagnated. Aborting." << std::endl;
			MMSP::Abort(-1);
		} else {
			if (rank==0)
				std::cout << "Converged in " << iter << " iterations." << std::endl;
		}

		output(initGrid,filename);
	}
}

template <int dim, typename T>
void update(grid<dim,vector<T> >& oldGrid, int steps)
{
	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

	grid<dim,vector<T> > newGrid(oldGrid);
	grid<dim,T> lapGrid(oldGrid, 0); // scalar grid from vector grid: zero fields

	// Make sure the grid spacing is correct
	for (int d=0; d<dim; d++) {
		dx(oldGrid,d) = deltaX;
		dx(newGrid,d) = deltaX;
		dx(lapGrid,d) = deltaX;

		// Set Neumann (zero-flux) boundary conditions
		if (x0(oldGrid, d) == g0(oldGrid, d)) {
			b0(oldGrid, d) = Neumann;
			b0(newGrid, d) = Neumann;
			b0(lapGrid, d) = Neumann;
		}
		if (x1(oldGrid, d) == g1(oldGrid, d)) {
			b1(oldGrid, d) = Neumann;
			b1(newGrid, d) = Neumann;
			b1(lapGrid, d) = Neumann;
		}
	}

	for (int step=0; step<steps; step++) {
		if (rank==0)
			print_progress(step, steps);

		ghostswap(oldGrid);

		for (int n=0; n<nodes(oldGrid); n++) {
			vector<int> x = position(oldGrid, n);
			const T& c = oldGrid(n)[0];
			const T& p = oldGrid(n)[1];
			lapGrid(n) = dfchemdc(c) + 2.0 * dfelecdc(p) - kappa*laplacian(oldGrid, x, 0);
		}

		ghostswap(lapGrid);

		// Apply Equation of Motion for composition (forward Euler)
		for (int n=0; n<nodes(oldGrid); n++) {
			vector<int> x = position(oldGrid, n);
			newGrid(n)[0] = oldGrid(n)[0] + dt*M*laplacian(lapGrid, x);
		}

		// Solve Poisson Equation for electrostatic potential
		unsigned int iter = RedBlackGaussSeidel(oldGrid, newGrid);
		if (iter >= max_iter) {
			if (rank==0)
				std::cerr << "Solver stagnated on step " << step << ". Aborting." << std::endl;
			MMSP::Abort(-1);
		}

		// ~ fin ~
		swap(oldGrid,newGrid);
	}
	ghostswap(oldGrid);
}

} // MMSP
#endif

#include"MMSP.main.hpp"
