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
	double fgrad = 0.0;
	double felec = 0.0;

	for (int n=0; n<nodes(GRID); n++) {
		vector<int> x = position(GRID, n);
		vector<double> gradc = gradient(GRID, x, 0);

		fchem += chemenergy(GRID(n)[0]);
		fgrad += gradc*gradc;
		felec += elecenergy(GRID(n)[0], GRID(n)[1]);
	}

	double F = dV*(fchem + 0.5*kappa*fgrad + felec); // Equation 5

	#ifdef MPI_VERSION
	double myF(F);
	MPI::COMM_WORLD.Allreduce(&myF, &F, 1, MPI_DOUBLE, MPI_SUM);
	#endif

	return F;
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
		GRID2D initGrid(2,0,200,0,200);
		for (int d=0; d<dim; d++)
			dx(initGrid,d) = deltaX;

		for (int n=0; n<nodes(initGrid); n++) {
			vector<int> x = position(initGrid, n);
			initGrid(n)[0] = cheminit(dx(initGrid,0)*x[0], dx(initGrid,1)*x[1]);
			initGrid(n)[1] = 0.0;
		}

		output(initGrid,filename);
	}
}

template <int dim, typename T>
void update(grid<dim,vector<T> >& oldGrid, int steps)
{
	grid<dim,vector<T> > newGrid(oldGrid);
	grid<dim,T> lapGrid(oldGrid, 0); // scalar grid from vector grid: zero fields
	// Make sure the grid spacing is correct
	for (int d=0; d<dim; d++) {
		dx(oldGrid,d) = deltaX;
		dx(newGrid,d) = deltaX;
		dx(lapGrid,d) = deltaX;
	}

	for (int step=0; step<steps; step++) {
		ghostswap(oldGrid);

		for (int n=0; n<nodes(oldGrid); n++) {
			vector<int> x = position(oldGrid, n);
			const T& c = oldGrid(n)[0];
			const T& p = oldGrid(n)[1];
			lapGrid(n) = dfdc(c) - kappa*laplacian(oldGrid, x, 0) + k*p;
		}

		ghostswap(lapGrid);

		// Apply Equation of Motion for composition (forward Euler)
		for (int n=0; n<nodes(oldGrid); n++) {
			vector<int> x = position(oldGrid, n);
			newGrid(n)[0] = oldGrid(n)[0] + dt*M*laplacian(lapGrid, x);
		}

		// Solve Poisson Equation for electrostatic potential (Gauss-Seidel)


		// ~ fin ~
		swap(oldGrid,newGrid);
	}
	ghostswap(oldGrid);
}

} // MMSP
#endif

#include"main.cpp"
