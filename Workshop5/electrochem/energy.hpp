// energy.hpp
// Energy functions for 2D and 3D Cahn-Hilliard model
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef CAHNHILLIARD_ENERGY
#define CAHNHILLIARD_ENERGY
#include<cmath>

// Composition parameters
const double Ca = 0.30; // alpha composition
const double Cb = 0.70; // beta composition
const double C0 = 0.50; // system composition
const double C1 = 0.04; // fluctuation magnitude

// Electrostatic parameters
const double k = 0.09;
const double epsilon = 9.0; // permittivity

// Physical parameters
const double rho = 5.0;
const double kappa = 2.0;
const double M = 5.0;

// Numerical parameters
const double CFL = 0.20;
const double deltaX = 1.0;
const double dt = std::pow(deltaX, 4)*CFL/(24.0*M*kappa);

double chemenergy(const double& C)
{
	// Equation 6
	const double A = C-Ca;
	const double B = Cb-C;
	return rho * A*A * B*B;
}

double elecenergy(const double& C, const double& P)
{
	// Equation 7
	return 0.5 * k * C * P;
}

double dfdc(const double& C)
{
	// d(chemenergy)/dc
	const double A = C-Ca;
	const double B = Cb-C;
	return 2.0 * rho * A * B * (Ca + Cb - 2.0 * C);
}

double cheminit(const double& x, const double& y)
{
	// Equation 12
	return C0 + C1 * ( std::cos(0.200*x)          * std::cos(0.11*y)
	                 + std::pow(std::cos(0.13*x)  * std::cos(0.087*y), 2.0)
	                 + std::cos(0.025*x - 0.15*y) * std::cos(0.07*x - 0.02*y));
}

#endif
