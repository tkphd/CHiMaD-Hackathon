// cahn-hilliard.hpp
// Algorithms for 2D and 3D Cahn-Hilliard model
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef CAHNHILLIARD_UPDATE
#define CAHNHILLIARD_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"cahn-hilliard.hpp"

const double q[2] = {0.1*std::sqrt(2.0), 0.1*std::sqrt(3.0)};
const double deltaX = 1.0;
const double Ca = 0.05;
const double Cb = 0.95;
const double Cm = 0.5*(Ca + Cb);
const double A = 2.0;
const double B = A/((Ca-Cm)*(Ca-Cm));
const double D = 2.0/(Cb-Ca);
const double K = 2.0;
const double dt = 0.005;
const double CFL = dt*K/std::pow(deltaX, 4);

double energydensity(double c)
{
	return -0.5*A*pow(c-Cm,2) + 0.25*B*pow(c-Cm,4) + 0.25*Ca*pow(c-Ca,4) + 0.25*Cb*pow(c-Cb,4);
}

namespace MMSP {

bool isOutside(const MMSP::vector<int>& x)
{
	if ((x[1] < 99) && (x[0]<40))
		return true;
	else if ((x[1] < 99) && (x[0]>59))
		return true;
	return false;
}

bool isBorderline(const MMSP::vector<int>& x)
{
	if ((x[1]==99) && (x[0]<40))
		return true;
	else if ((x[1]==99) && (x[0]>59))
		return true;
	else if ((x[1]<99) && (x[0]==40))
		return true;
	else if ((x[1]<99) && (x[0]==60))
		return true;
	return false;
}

// custom Laplacian for boundary points
template <int dim, typename T>
T zfLaplacian(const grid<dim, T>& GRID, const vector<int>& x)
{
  T laplacian = 0.0;
  MMSP::vector<int> s = x;
  const T& y = GRID(x);

  for (int i=0; i<dim; i++) {
    s[i] += 1;
    const T& yh = GRID(s);
    s[i] -= 2;
    const T& yl = GRID(s);
    s[i] += 1;

    double weight = 1.0 / (dx(GRID, i) * dx(GRID, i));
    if (x[1]==99 && x[0]==40) // low side
    	laplacian += weight * (yh - y);
    else if (x[1]==99 && x[0]==60) // high side
    	laplacian += weight * (-y + yl);
    else
    	laplacian += weight * (yh - 2.0 * y + yl);
  }
  return laplacian;
}


void generate(int dim, const char* filename)
{
	if (dim!=2) {
		std::cerr<<"ERROR: CHiMaD problems are 2-D, only!"<<std::endl;
		std::exit(-1);
	}

	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

	if (dim==2) {
		MMSP::grid<2,double> grid(1,0,100,0,120);
		for (int d=0; d<dim; d++){
			dx(grid,d) = deltaX;
			if (MMSP::x0(grid,d)==MMSP::g0(grid,d))
			    MMSP::b0(grid,d) = Neumann; // enumerated in MMSP.utility.hpp
			else if (MMSP::x1(grid,d)==MMSP::g1(grid,d))
			    MMSP::b1(grid,d) = Neumann; // enumerated in MMSP.utility.hpp
		}

		for (int i=0; i<nodes(grid); i++) {
			MMSP::vector<int> x = position(grid,i);
			if (isOutside(x))
				grid(x) = 0.0;
			else
				grid(x) = 0.45 + 0.01 * std::cos(x[0]*dx(grid,0)*q[0] + x[1]*dx(grid,1)*q[1]);
		}

		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
		output(grid,filename);
		if (rank==0)
			std::cout<<"Timestep is "<<dt<<" (CFL="<<CFL<<')'<<std::endl;
	}
}

template <int dim, typename T>
void update(MMSP::grid<dim,T>& grid, int steps)
{
	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

	// Make sure the grid spacing is correct
	for (int d=0; d<dim; d++) {
		dx(grid,d) = deltaX;
		if (MMSP::x0(grid,d)==MMSP::g0(grid,d))
		    MMSP::b0(grid,d) = Neumann; // enumerated in MMSP.utility.hpp
		else if (MMSP::x1(grid,d)==MMSP::g1(grid,d))
		    MMSP::b1(grid,d) = Neumann; // enumerated in MMSP.utility.hpp
  }

    // Let's be absolutely explicit about BCs here.
	MMSP::grid<dim,T> update(grid);
	for (int d=0; d<dim; d++) {
		dx(update,d) = deltaX;
		if (MMSP::x0(update,d)==MMSP::g0(update,d))
		    MMSP::b0(update,d) = Neumann; // enumerated in MMSP.utility.hpp
		else if (MMSP::x1(update,d)==MMSP::g1(update,d))
		    MMSP::b1(update,d) = Neumann; // enumerated in MMSP.utility.hpp
  }
	MMSP::grid<dim,T> temp(grid);
	for (int d=0; d<dim; d++) {
		dx(temp,d) = deltaX;
		if (MMSP::x0(temp,d)==MMSP::g0(temp,d))
		    MMSP::b0(temp,d) = Neumann; // enumerated in MMSP.utility.hpp
		else if (MMSP::x1(temp,d)==MMSP::g1(temp,d))
		    MMSP::b1(temp,d) = Neumann; // enumerated in MMSP.utility.hpp
  }


	for (int step=0; step<steps; step++) {
		for (int n=0; n<nodes(grid); n++) {
			MMSP::vector<int> x = position(grid,n);
			if (isOutside(x))
				continue;
			double c = grid(x);
			double dfdc = -A*(c-Cm) + B*pow(c-Cm, 3) + Ca*pow(c-Ca, 3) + Cb*pow(c-Cb, 3);
			temp(x) = dfdc - K*zfLaplacian(grid,x);
		}
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
		ghostswap(temp);

		double energy = 0.0;
		for (int n=0; n<nodes(grid); n++) {
			MMSP::vector<int> x = position(grid,n);
			if (isOutside(x))
				continue;
			update(x) = grid(x)+dt*D*zfLaplacian(temp,x);
			energy += dx(grid)*dy(grid)*energydensity(update(x));
		}
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		double myenergy = energy;
		MPI::COMM_WORLD.Allreduce(&myenergy, &energy, 1, MPI_DOUBLE, MPI_SUM);
		#endif
		if (rank==0)
			std::cout<<energy<<std::endl;

		swap(grid,update);
		ghostswap(grid);
	}
}

} // MMSP
#endif

#include"MMSP.main.hpp"
