/** stokes.cpp
 ** Algorithms for 2D Stokes model with convex splitting.
 **
 ** This is an advanced example, with precompiler directives switching between
 ** shared-memory (OpenMP) and message-passing (MPI) parallelism, with the equations
 ** of motion discretized using semi-implicit convex splitting for guaranteed energy
 ** stability, i.e. dE/dt<=0, and with successive over-relaxation to accelerate
 ** convergence.
 **
 ** Questions/comments to trevor.keller@gmail.com (Trevor Keller)
 **/

#ifndef STOKES_UPDATE
#define STOKES_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include<cfloat>
#include"stokes.hpp"

// Model parameters
const double rho = 100.;
const double mu = 1.;
vector<double> g(2, 0.); g[1] = -0.001;

inline double profile(double y)
{
	return 0.009 - 0.001 * std::pow(y - 3.0, 2);
}

// Numerical constants
double tolerance = 1.0e-12;		// Choose wisely. 1e-10 is the minimum toloerance for which mass is conserved.
								// Tighter tolerance is better, but increases runtime.
unsigned int residual_step = 5;	// number of iterations between residual computations

template<typename T> inline
double energy_density(const T& C)
{
	// Minima at C = +/- 1, local maximum at C=0
	return 0.25*pow(C,4.0) - 0.5*pow(C,2.0);
}
template<typename T> inline
double full_dfdc(const T& C)
{
	return pow(C,3.0) - C;
}
template<typename T> inline
double contractive_dfdc(const T& C)
{
	return pow(C,3.0);
}
template<typename T> inline
double nonlinear_coeff(const T& C)
{
	return pow(C,2.0);
}
template<typename T> inline
double expansive_dfdc(const T& C)
{
	return -C;
}

namespace MMSP
{

// Discrete Laplacian operator missing the central value, for implicit source terms
template<int dim, typename T>
double fringe_laplacian(const grid<dim,vector<T> >& GRID, const vector<int>& x, const int field)
{
	double laplacian = 0.0;
	vector<int> s = x;

	for (int i=0; i<dim; i++) {
		s[i] += 1;
		const T& yh = GRID(s)[field];
		s[i] -= 2;
		const T& yl = GRID(s)[field];
		s[i] += 1;

		double weight = 1.0 / pow(dx(GRID, i),2.0);
		laplacian += weight * (yh + yl);
	}
	return laplacian;
}

void generate(int dim, const char* filename)
{
	if (dim==2) {
		// Mesh size
		double Lx = 30;
		double Ly = 6;

		// Mesh resolution
		double h = 1.0;

		GRID2D initGrid(3,0,std::floor(Lx/h),0,std::floor(Ly/h)); // field 0 is c, field 1 is mu
		for (int d=0; d<dim; d++)
			dx(initGrid,d) = h;

		vector<double> blank(dim+1, 0.0);
		bool touchedCorner = false;
		for (int n=0; n<nodes(initGrid); n++) {
			vector<int> x = position(initGrid,n);
			initGrid(n) = blank;

			// Set pressure: interpolate between inlet at 1 and outlet at zero, so grad(p) is not zero
			initGrid(n)[dim] = 1.0 - double(x[0] - g0(initGrid,0))/(g1(initGrid,0) - g0(initGrid,0));

			if (x[0] < g1(initGrid, 0)-1 && x[1] > g0(initGrid, 1) && x[1] < g1(initGrid, 1)-1) {
				// point is within the domain
				double v0 = profile(dx(initGrid,1) * x[1]);
				// interpolate between inlet at v0 and outlet at, initially, zero
				// (other choices of initial condition may be more reasonable)
				initGrid(n)[0] = v0 * (1.0 - double(x[0] - g0(initGrid,0))/(g1(initGrid,0) - g0(initGrid,0)));

			}
		}

		for (int n=0; n<nodes(initGrid); n++) {
			initGrid(n)[1] = full_dfdc(initGrid(n)[0]) - K * laplacian(initGrid, x, 0);
		}

		output(initGrid,filename);
		//reportEnergy(initGrid);
	} else {
		std::cerr << "Error: " << dim << "-D code has not been implemented."<<std::endl;
		MMSP::Abort(-1);
	}
}

template<int dim, typename T>
bool update(grid<dim,vector<T> >& oldGrid, int steps)
{
	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

	ghostswap(oldGrid);

	grid<dim,vector<T> > newGrid(oldGrid);   // new values at each point and initial guess for iteration

	newGrid.copy(oldGrid); // deep copy: includes data and ghost cells. Expensive.

	// Make sure the grid spacing is correct.
	for (int d=0; d<dim; d++) {
		dx(oldGrid,d) = deltaX;
		dx(newGrid,d) = deltaX;
	}

	double gridSize = static_cast<double>(nodes(oldGrid));
	#ifdef MPI_VERSION
	double localGridSize = gridSize;
	MPI::COMM_WORLD.Allreduce(&localGridSize, &gridSize, 1, MPI_DOUBLE, MPI_SUM);
	#endif

	double lapWeight = 0.0;
	double dV = 1.0;
	for (int d=0; d<dim; d++) {
		lapWeight += 2.0/pow(dx(oldGrid,d),2.0); // dim=2 -> lapWeight = 4/h^2 if dy=dx=h
		dV *= dx(oldGrid,d);
	}

	static double residual = 1.0;
	bool evolving = true;

	for (int iter=0; iter<steps && evolving; iter++) {
		if (rank==0)
			print_progress(iter, steps);

		ghostswap(oldGrid);
		ghostswap(newGrid);

		/*  ==== RED-BLACK GAUSS SEIDEL ====
		    Iterate over a checkerboard, updating first red then black tiles.
		    This method eliminates the third "guess" grid, and should converge faster.
		    In 2-D and 3-D, if the sum of indices is even, then the tile is Red; else, Black.
		*/

		for (int color=1; color>-1; color--) {
			// If color==1, skip BLACK tiles, which have Σx[d] odd
			// If color==0, skip RED tiles, which have Σx[d] even

			// OpenMP parallel loop over nodes
			#pragma omp parallel for schedule(dynamic)
			for (int n=0; n<nodes(oldGrid); n++) {
				vector<int> x = position(oldGrid,n);
				int x_sum=0;
				for (int d=0; d<dim; d++)
					x_sum += x[d];
				if (x_sum%2 == color)
					continue;

				T cOld = oldGrid(n)[0];
				T cGuess = newGrid(n)[0]; // value from last "guess" iteration
				T uGuess = newGrid(n)[1];

				// A is defined by the last guess, stored in newGrid(n). It is a 2x2 matrix.
				//T A11 = 1.0;
				T A12 = dt*D*lapWeight;
				T A21 = -nonlinear_coeff(newGrid(n)[0]) - K*lapWeight;
				//T A22 = 1.0;

				T detA = 1.0 - (A12*A21); // determinant of A

				// B is defined by the last value, stored in oldGrid(n), and the last guess, stored in newGrid(n). It is a 2x1 column.
				T B1 = cOld + D*dt*fringe_laplacian(newGrid, x, 1);
				T B2 = expansive_dfdc(cOld) - K*fringe_laplacian(newGrid, x, 0);

				// Solve the iteration system AX=B using Cramer's rule
				T cNew = (B1 - B2*A12)/detA; // X1
				T uNew = (B2 - B1*A21)/detA; // X2

				// (Don't) Apply relaxation
				newGrid(n)[0] = omega*cNew + (1.0 - omega)*cGuess;
				newGrid(n)[1] = omega*uNew + (1.0 - omega)*uGuess;

			}
			ghostswap(newGrid);   // fill in the ghost cells; does nothing in serial


			evolving = (residual>tolerance);

			iter++;

			/*  ==== RESIDUAL ====
			    The residual is computed from the original matrix form, Ax=b:
			    any Old term goes into B, while any New term goes in AX. Note that
			    this is not the iteration matrix, it is the original system of equations.
			*/

			if (iter<residual_step || iter%residual_step==0) {
				double normB = 0.0;
				residual = 0.0;

				// OpenMP parallel loop over nodes
				#pragma omp parallel for schedule(dynamic)
				for (int n=0; n<nodes(oldGrid); n++) {
					vector<int> x = position(oldGrid,n);
					T lapC = laplacian(newGrid, x, 0);
					T lapU = laplacian(newGrid, x, 1);

					T cOld = oldGrid(n)[0];
					T cNew = newGrid(n)[0];
					T uNew = newGrid(n)[1];

					T B1 = cOld;
					T B2 = expansive_dfdc(cOld);

					T AX1 = cNew - D*dt*lapU;
					T AX2 = uNew - contractive_dfdc(cNew) + K*lapC;

					// Compute the Error from parts of the solution
					double R1 = B1 - AX1;
					double R2 = B2 - AX2;

					double error = R1*R1 + R2*R2;
					#pragma omp critical
					{
						residual += error;
						normB += B1*B1 + B2*B2;
					}
				}

				#ifdef MPI_VERSION
				double localResidual=residual;
				MPI::COMM_WORLD.Allreduce(&localResidual, &residual, 1, MPI_DOUBLE, MPI_SUM);
				double localNormB=normB;
				MPI::COMM_WORLD.Allreduce(&localNormB, &normB, 1, MPI_DOUBLE, MPI_SUM);
				#endif

				residual = sqrt(residual/normB)/(2.0*gridSize);
			}

		}

		swap(oldGrid, newGrid);
	}
	if (rank==0)
		std::cout<<std::flush;

	return evolving;
}

} // MMSP
#endif

#include"main.cpp"
