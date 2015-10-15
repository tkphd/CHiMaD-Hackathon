// ostwald.hpp
// Ostwald ripening algorithms for 2D and 3D phase field methods
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef OSTWALD_UPDATE
#define OSTWALD_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"ostwald.hpp"

const double deltaX = 1.0;
const double Ca = 0.05;
const double Cb = 0.95;
const double Cm = 0.5*(Ca+Cb);
const double A = 2.0;
const double B = A/std::pow(Ca-Cm,2);
const double g = 2.0/std::pow(Cb-Ca,2);
const double delta = 1.0;
const double epsilon = 3.0;
const double Da = g/std::pow(delta,2);
const double Db = g/std::pow(delta,2);
const double kappa = 2.0;
const double D = 1.0;
const double L = 1.0;

const double dt = 0.005; //std::std::pow(deltaX, 4.0)/(160.0 * K); // might need updating
const double epsi[11] = {0.0, 0.979285, 0.219812, 0.837709, 0.695603, 0.225115, 0.389266, 0.585953, 0.614471, 0.918038, 0.518569};
const double CFL = dt*D*kappa/std::pow(deltaX,4);

double energydensity(const MMSP::vector<double>& value)
{
  double C = value[0];
  double f1 = -0.5*A*pow(C-Cm,2) + 0.25*B*pow(C-Cm,4) + 0.25*Da*pow(C-Ca,4) + 0.25*Db*pow(C-Cb,4);
  double f2 = 0.0;
  for (int i=1; i<length(value); i++)
    f2+=-0.5*g*pow(C-Ca,2)*pow(value[i],2)+0.25*delta*pow(value[i],4);
  double f3 = 0.0;
  for (int i=1; i<length(value); i++)
    for (int j=i+1; j<length(value); j++)
      f3 += epsilon*pow(value[i],2)*pow(value[j],2); // Half vanishes to enable double-counting
  return f1 + f2 + f3;
}

namespace MMSP{

bool isOutside(const MMSP::vector<int>& x)
{
  if ((x[1]<99) && ((x[0]<40) || (x[0]>59)))
    return true;
  return false;
}

bool isBorderline(const MMSP::vector<int>& x)
{
  if ((x[1]==99) && (x[0]>40) && (x[0]<60))
    return true;
  else if ((x[1]<99) && ((x[0]==40) || (x[0]==59)))
    return true;

  return false;
}

// custom Laplacian for boundary points and specified field
template <int dim, typename T>
T zfLaplacian(const grid<dim, MMSP::vector<T> >& GRID, const vector<int>& x, const int field)
{
  T laplacian = 0.0;
  MMSP::vector<int> s = x;
  const T& y = GRID(x)[field];

  for (int i=0; i<dim; i++) {
    s[i] += 1;
    const T& yh = GRID(s)[field];
    s[i] -= 2;
    const T& yl = GRID(s)[field];
    s[i] += 1;

    double weight = 1.0 / (dx(GRID, i) * dx(GRID, i));
    if (i==0 && x[1]<99 && x[0]==40) // low side vacant
        laplacian += weight * (yh - y);
    else if (i==0 && x[1]<=99 && x[0]==59) // high side vacant
        laplacian += weight * (-y + yl);
    else if (i==1 && x[1]==99 && (x[0]<=40 || x[0]>=59))
        laplacian += weight * (yh - y);
    else
      laplacian += weight * (yh - 2.0 * y + yl);
  }
  return laplacian;
}

template <int dim, typename T>
T zfLaplacian(const grid<dim, MMSP::scalar<T> >& GRID, const vector<int>& x)
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
    if (i==0 && x[1]<99 && x[0]==40) // low side vacant
        laplacian += weight * (yh - y);
    else if (i==0 && x[1]<=99 && x[0]==59) // high side vacant
        laplacian += weight * (-y + yl);
    else if (i==1 && x[1]==99 && (x[0]<=40 || x[0]>=59))
        laplacian += weight * (yh - y);
    else
      laplacian += weight * (yh - 2.0 * y + yl);
  }
  return laplacian;
}

void generate(int dim, const char* filename)
{
	if (dim!=2) {
		std::cerr<<"ERROR: CHiMaD benchmarks are 2-D, only!"<<std::endl;
		std::exit(-1);
	}

	int rank=0;
  #ifdef MPI_VERSION
  rank = MPI::COMM_WORLD.Get_rank();
  #endif

	const double q[2] = {0.1*std::sqrt(2.0), 0.1*std::sqrt(3.0)};
	double qi[11][2];
	qi[0][0] = 0.0;
	qi[0][1] = 0.0;
	for (int i=1; i<11; i++) {
		qi[i][0] = 0.01*std::sqrt(23+i);
		qi[i][1] = 0.01*std::sqrt(149+i);
	}
	if (dim==2) {
		MMSP::grid<2,MMSP::vector<double> > grid(11,0,100,0,120);

    for (int d=0; d<dim; d++) {
      dx(grid,d) = deltaX;
			if (MMSP::x0(grid,d)==MMSP::g0(grid,d))
          MMSP::b0(grid,d) = Neumann; // enumerated in MMSP.utility.hpp
      else if (MMSP::x1(grid,d)==MMSP::g1(grid,d))
          MMSP::b1(grid,d) = Neumann; // enumerated in MMSP.utility.hpp
		}

    for (int i=0; i<nodes(grid); i++) {
      MMSP::vector<int> x = position(grid,i);
      if (isOutside(x)) {
      	grid(x)[0] = 0.0;
				for (int i=1; i<fields(grid); i++)
		      grid(x)[i] = 0.0;
      } else {
      	grid(x)[0] = 0.5 + 0.01 * std::cos(x[0]*dx(grid,0)*q[0] + x[1]*dx(grid,1)*q[1]); // conc
				for (int i=1; i<fields(grid); i++)
	      	grid(x)[i] = 0.0 + 0.01 * epsi[i] * std::pow(std::cos(x[0]*dx(grid,0)*qi[i][0] + x[1]*dx(grid,1)*qi[i][1]),2); // phase
	    }
    }

    #ifdef MPI_VERSION
    MPI::COMM_WORLD.Barrier();
    #endif
		MMSP::output(grid,filename);
		if (rank==0)
      std::cout<<"Timestep is "<<dt<<" (CFL="<<CFL<<')'<<std::endl;
	}
}

void update(MMSP::grid<1,MMSP::vector<double> >& grid, int steps)
{
	std::cerr<<"ERROR: CHiMaD benchmarks are 2-D, only!"<<std::endl;
	std::exit(-1);
}

void update(MMSP::grid<3,MMSP::vector<double> >& grid, int steps)
{
	std::cerr<<"ERROR: CHiMaD benchmarks are 2-D, only!"<<std::endl;
	std::exit(-1);
}

void update(MMSP::grid<2,MMSP::vector<double> >& grid, int steps)
{
	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

  for (int d=0; d<2; d++) {
    dx(grid,d) = deltaX;
		if (MMSP::x0(grid,d)==MMSP::g0(grid,d))
        MMSP::b0(grid,d) = Neumann; // enumerated in MMSP.utility.hpp
    else if (MMSP::x1(grid,d)==MMSP::g1(grid,d))
        MMSP::b1(grid,d) = Neumann; // enumerated in MMSP.utility.hpp
	}
	ghostswap(grid);

	MMSP::grid<2,MMSP::vector<double> > update(grid);
  for (int d=0; d<2; d++) {
    dx(update,d) = deltaX;
		if (MMSP::x0(update,d)==MMSP::g0(update,d))
        MMSP::b0(update,d) = Neumann; // enumerated in MMSP.utility.hpp
    else if (MMSP::x1(update,d)==MMSP::g1(update,d))
        MMSP::b1(update,d) = Neumann; // enumerated in MMSP.utility.hpp
	}

	MMSP::grid<2,MMSP::scalar<double> > wspace(grid,1);
  for (int d=0; d<2; d++) {
    dx(wspace,d) = deltaX;
		if (MMSP::x0(wspace,d)==MMSP::g0(wspace,d))
        MMSP::b0(wspace,d) = Neumann; // enumerated in MMSP.utility.hpp
    else if (MMSP::x1(wspace,d)==MMSP::g1(wspace,d))
        MMSP::b1(wspace,d) = Neumann; // enumerated in MMSP.utility.hpp
	}


	for (int step=0; step<steps; step++) {
		for (int n=0; n<nodes(grid); n++) {
			MMSP::vector<int> x=position(grid,n);
			double sum = 0.0;
			for (int i=1; i<fields(grid); i++)
				sum += std::pow(grid(x)[i],2);

			double C = grid(x)[0];
			double lap = zfLaplacian(grid, x, 0);

			wspace(x) = -A*(C-Cm) + B*std::pow(C-Cm,3) + Da*std::pow(C-Ca,3) + Db*std::pow(C-Cb,3) - g*(C-Ca)*sum - kappa*lap;
		}
		ghostswap(wspace);

		double energy = 0.0;
		int err=0;
		for (int n=0; n<nodes(grid); n++) {
			MMSP::vector<int> x=position(grid,n);
			if (1) {
				double lap = zfLaplacian(wspace, x);
				update(x)[0] = grid(x)[0]+dt*D*lap;
			}

			double C = grid(x)[0];
			double sum = 0.0;
			for (int i=1; i<fields(grid); i++)
				sum += std::pow(grid(x)[i],2);

			for (int i=1; i<fields(grid); i++) {
				double phase = grid(x)[i];
				double lap = zfLaplacian(grid, x, i);

				update(x)[i] = grid(x)[i] + dt*L*(g*std::pow(C-Ca,2)*phase - delta*std::pow(phase,3)
				                                      - epsilon*phase*(sum-std::pow(phase,2)) + kappa*lap);
				double de = energydensity(update(x));
				if (std::isfinite(de))
					energy += dx(grid)*dy(grid)*energydensity(update(x));
				else
					err++;
			}
		}
		#ifdef MPI_VERSION
		double myEnergy=energy;
		double myErr=err;
		MPI::COMM_WORLD.Allreduce(&myEnergy,&energy,1,MPI_DOUBLE,MPI_SUM);
		MPI::COMM_WORLD.Allreduce(&myErr,&err,1,MPI_DOUBLE,MPI_SUM);
		#endif
		if (rank==0)
			std::cout<<energy<<std::endl;
		if (rank==0 && err>0)
			std::cerr<<err<<" invalid pixels"<<std::endl;
		swap(grid,update);
		ghostswap(grid);
	}
}

double zfLaplacian(const grid<2, vector<double> >& GRID, const vector<int>& x, const int field)
{
  double laplacian(0.0);
  vector<int> s = x;

  const double& y = GRID(x)[field];

  for (int i=0; i<2; i++) {
    s[i] += 1;
    const double& yh = GRID(s)[field];
    s[i] -= 2;
    const double& yl = GRID(s)[field];
    s[i] += 1;

    double weight = 1.0 / (dx(GRID, i) * dx(GRID, i));
    laplacian += weight * (yh - 2.0 * y + yl);
  }
  return laplacian;
}


} // namespace MMSP

#endif

#include"MMSP.main.hpp"
