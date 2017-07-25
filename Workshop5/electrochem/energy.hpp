// energy.hpp
// Energy functions for 2D and 3D Cahn-Hilliard model
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef CAHNHILLIARD_ENERGY
#define CAHNHILLIARD_ENERGY
#include<cmath>

#define cid 0
#define uid 1
#define pid 2

// Composition parameters
const double Ca = 0.30; // alpha composition
const double Cb = 0.70; // beta composition
const double C0 = 0.50; // system composition
const double C1 = 0.04; // fluctuation magnitude

// Electrostatic parameters
const double k = 0.09;
const double epsilon = 90.0; // permittivity

// Physical parameters
const double rho = 5.0;
const double kappa = 2.0;
const double M = 5.0;

// Numerical parameters
const double CFL = 0.25;
const double deltaX = 1.0;
const double dt = std::pow(deltaX, 4)*CFL/(24.0*M*kappa);

// Gauss-Seidel parameters
double tolerance = 1.0e-12;		// Choose wisely. 1e-10 is the minimum toloerance for which mass is conserved.
unsigned int residual_step = 5;	// number of iterations between residual computations
unsigned int max_iter = 10000;	// don't let the solver stagnate
double omega = 1.20;			// relaxation parameter for SOR. omega=1 is stock Gauss-Seidel.


// Energy equations
double cheminit(const double& x, const double& y)
{
	// Equation 12
	return C0 + C1 * ( std::cos(0.200*x)          * std::cos(0.11*y)
	                 + std::pow(std::cos(0.13*x)  * std::cos(0.087*y), 2.0)
	                 + std::cos(0.025*x - 0.15*y) * std::cos(0.07*x - 0.02*y));
}

template<typename T>
double chemenergy(const T& C)
{
	// Equation 6
	const double A = C-Ca;
	const double B = Cb-C;
	return rho * A*A * B*B;
}

template<typename T>
double elecenergy(const T& C, const T& P)
{
	// Equation 7
	return 0.5 * k * C * P;
}


// Energy derivatives
template<typename T>
double dfchemdc(const T& C)
{
	// d(chemenergy)/dc
	const double A = C-Ca;
	const double B = Cb-C;
	return 2.0 * rho * A * B * (Ca + Cb - 2.0 * C);
}

template<typename T>
double dfelecdc(const T& P)
{
	return 0.5 * k * P;
}

template<typename T>
double dfcontractivedc(const T& C, const T& Cnew)
{
	return 2.0 * rho * Cnew * (2.0 * C*C + Ca*Ca + 4.0*Ca*Cb + Cb*Cb) ;
}

template<typename T>
double dfexpansivedc(const T& C, const T& P)
{
	return -2.0 * rho * (3.0 * C*C * (Ca + Cb) - Ca * Cb * (Ca + Cb)) + 2.0 * dfelecdc(P);
}

// Discrete Laplacian operator missing the central value, for implicit source terms
template<int dim, typename T>
double fringe_laplacian(const MMSP::grid<dim,MMSP::vector<T> >& GRID, const MMSP::vector<int>& x, const int field)
{
	double laplacian = 0.0;
	MMSP::vector<int> s = x;

	for (int i=0; i<dim; i++) {
		s[i] += 1;
		const T& yh = GRID(s)[field];
		s[i] -= 2;
		const T& yl = GRID(s)[field];
		s[i] += 1;

		double weight = 1.0 / pow(dx(GRID, i), 2.0);
		laplacian += weight * (yh + yl);
	}
	return laplacian;
}

#endif
