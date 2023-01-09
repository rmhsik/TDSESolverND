#include <cmath>
#include "potential.h"


cdouble potential( double i, double k, double t, Hamiltonian *ham){
}

//cdouble potential(double i, double k, double t, Hamiltonian *ham){
//    double bfield_i = ham->_param->B0i;
//    return -1.0/sqrt(i*i+k*k + 0.65) + 0.125*bfield_i*bfield_i*(i*i+k*k);
//}

cdouble potential_XZ( double i, double k, double ti, Hamiltonian *ham){
    return -1.0/sqrt(i*i + k*k + 0.65);
}

