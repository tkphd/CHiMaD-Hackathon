// energy.hpp
// Ostwald ripening algorithms for 2D and 3D phase field methods
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef OSTWALD_ENERGY
#define OSTWALD_ENERGY
#include<cmath>

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
const double CFL = 0.125;
const double dt = std::pow(deltaX,4)*CFL/(32.0*D*kappa);

double energydensity(const MMSP::vector<double>& value)
{
	const double& C = value[0];
	double f = -0.5*A*std::pow(C-Cm,2) + 0.25*B*std::pow(C-Cm,4) + 0.25*Da*std::pow(C-Ca,4) + 0.25*Db*std::pow(C-Cb,4);
	double sum = 0.0;

	for (int i=1; i<length(value); i++)
	    sum += std::pow(value[i],2);

	for (int i=1; i<length(value); i++) {
	    double psq = std::pow(value[i],2);
   		f +=-0.5*g*std::pow(C-Ca,2)*psq + 0.25*delta*psq*psq // f2
		   + 0.5*epsilon*psq*(sum - psq); // f3
	}

	return f;
}

double df1dc(const double& C)
{
 return -A*(C-Cm) + B*std::pow(C-Cm,3) + Da*std::pow(C-Ca,3) + Db*std::pow(C-Cb,3);
}

double df2dc(const double& C, const double& sum)
{
 return -g*(C-Ca)*sum;
}

double df2deta(const double& C, const double& phase)
{
    return -1.0*g*std::pow(C-Ca,2)*phase + delta*std::pow(phase,3);
}

double df3deta(const double& phase, const double& sum)
{
    return epsilon*phase*(sum-std::pow(phase,2));
}

#endif
