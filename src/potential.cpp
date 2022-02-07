#include <cmath>
#include "potential.h"


cdouble potential( double i, double k, double t, Hamiltonian *ham){
    cdouble omegaB = ham->_param->B0k/(sqrt(8));
    cdouble omega0 = ham->_param->w0;
    cdouble phiBk = ham->_param->phiBk;
    return 0.5*omegaB*omegaB*i*i - 0.5*omegaB*omegaB*i*i*cos(2.0*omega0*t + phiBk);
}

//cdouble potential(double i, double k, double t, Hamiltonian *ham){
//    double bfield_i = ham->_param->B0i;
//    return -1.0/sqrt(i*i+k*k + 0.65) + 0.125*bfield_i*bfield_i*(i*i+k*k);
//}

cdouble potential_X( double i, double k, double ti, Hamiltonian *ham){
    return -1.0/sqrt(i*i + k*k + 2.0);
}
cdouble potential_XZ( double i, double k, double ti, Hamiltonian *ham){
    return -1.0/sqrt(i*i + k*k + 0.65);
}
cdouble potential_RZ( double i, double k, double ti, Hamiltonian *ham){
    return -1.0/sqrt(i*i + k*k );
}
