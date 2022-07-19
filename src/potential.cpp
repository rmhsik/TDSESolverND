#include <cmath>
#include "potential.h"


cdouble potential( double i, double j, double k, double t, Hamiltonian *ham){
    cdouble omegaB = ham->_param->B0k/(sqrt(8));
    cdouble omega0 = ham->_param->w0;
    cdouble phiBk = ham->_param->phiBk;
    cdouble phiEk = ham->_param->phiEk;
    double E0 = ham->_param->E0k;
    double ti = 1.8*ham->_param->period;
    return 0.5*omegaB*omegaB*i*i - 0.5*omegaB*omegaB*i*i*cos(2.0*omega0*(t+ti) + phiBk) - E0*sin(omega0*(t+ti)+phiEk)*i;
}

//cdouble potential(double i, double k, double t, Hamiltonian *ham){
//    double bfield_i = ham->_param->B0i;
//    return -1.0/sqrt(i*i+k*k + 0.65) + 0.125*bfield_i*bfield_i*(i*i+k*k);
//}

cdouble potential_X( double i, double j, double k, double ti, Hamiltonian *ham){
    return -1.0/sqrt(i*i + k*k + 2.0);
}
cdouble potential_XZ( double i, double j, double k, double ti, Hamiltonian *ham){
    return -1.0/sqrt(i*i + k*k + 0.65);
}
cdouble potential_RZ( double i, double j, double k, double ti, Hamiltonian *ham){
    return -1.0/sqrt(i*i + k*k );
}
cdouble potential_XYZ( double i, double j, double k, double ti, Hamiltonian *ham){
    return -1.0/sqrt(i*i + k*k + j*j);
}

cdouble potential_argon_XYZ (double i, double j, double k, double ti, Hamiltonian *ham){
    // doi: 10.1088/2399-6528/ab9a68
    double r  = sqrt(i*i + j*j + k*k + 0.4);
    double r_inv = 1.0/r;
    //const double c  =   0.8103;
    //const double c0 =   1.0000;  const double zc =  17.0000;
    //const double a1 = -15.9583;  const double b1 =   1.2305;
    //const double a2 = -27.7467;  const double b2 =   4.3946;
    //const double a3 =   2.1768;  const double b3 =  90.0460;
    //return -c0/r - zc*exp(-c*r)/r - a1*exp(-b1*r) - a2*exp(-b2*r) - a3*exp(-b3*r);
    //
    // Tong: https://doi.org/10.1088/0953-4075/38/15/001
    double zc =   1.0;     double a1 =  16.039;
    double a2 =   2.007;   double a3 = -25.543;
    double a4 =   4.525;   double a5 =   0.961;
    double a6 =   0.443;
    return -r_inv*(zc + a1*exp(-a2*r) + a3*r*exp(-a4*r) + a5*exp(-a6*r));

    //return -r_inv;//*(1+5.4*exp(-r)+(17-5.4)*exp(-3.682*r));

}
