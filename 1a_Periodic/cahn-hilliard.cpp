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
const double dt = 0.005; //std::pow(deltaX, 4)/(160 * K); // 0.003125 appears stable

double energydensity(double c)
{
	return -0.5*A*pow(c-Cm,2) + 0.25*B*pow(c-Cm,4) + 0.25*Ca*pow(c-Ca,4) + 0.25*Cb*pow(c-Cb,4);
}

namespace MMSP {

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
		MMSP::grid<2,double> grid(1,0,200,0,200);
		for (int d=0; d<dim; d++)
			dx(grid,d) = deltaX;

		for (int i=0; i<nodes(grid); i++) {
			MMSP::vector<int> x = position(grid,i);
			grid(x) = 0.45 + 0.01 * std::cos(x[0]*dx(grid,0)*q[0] + x[1]*dx(grid,1)*q[1]);
		}

		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
		output(grid,filename);
		if (rank==0)
			std::cout<<"Timestep is "<<dt<<" (Co~1/40)"<<std::endl;
	}
}

template <int dim, typename T>
void update(MMSP::grid<dim,T>& grid, int steps)
{
	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif
	MMSP::grid<dim,T> update(grid);
	MMSP::grid<dim,T> temp(grid);
	// Make sure the grid spacing is correct
	for (int d=0; d<dim; d++) {
		dx(grid,d) = deltaX;
		dx(update,d) = deltaX;
		dx(temp,d) = deltaX;
	}


	for (int step=0; step<steps; step++) {
		for (int i=0; i<nodes(grid); i++) {
			MMSP::vector<int> x = position(grid,i);
			double c = grid(x);
			double dfdc = -A*(c-Cm) + B*pow(c-Cm, 3) + Ca*pow(c-Ca, 3) + Cb*pow(c-Cb, 3);
			temp(i) = dfdc - K*laplacian(grid,x);
		}
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
		ghostswap(temp);

		double energy = 0.0;
		for (int i=0; i<nodes(grid); i++) {
			MMSP::vector<int> x = position(grid,i);
			update(x) = grid(x)+dt*D*laplacian(temp,x);
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
