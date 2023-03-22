#include <cmath>
#include "potential.h"


cdouble potential( double i, double k, double t, Hamiltonian *ham){
}

//cdouble potential(double i, double k, double t, Hamiltonian *ham){
//    double bfield_i = ham->_param->B0i;
//    return -1.0/sqrt(i*i+k*k + 0.65) + 0.125*bfield_i*bfield_i*(i*i+k*k);
//}

cdouble potential_hydrogen_XZ(double i, double k, double ti, Hamiltonian *ham){
    return -1.0/sqrt(i*i + k*k + 0.65);
}

cdouble potential_helium_XZ(double i, double k, double ti, Hamiltonian *ham){
    double r = sqrt(i*i +k*k);
    double inv_r = 1.0/sqrt(i*i + k*k+0.05);
    double a[7] = {1.0, 1.2231, 0.662, -1.325, 1.236, -0.231, 0.480};
    return -(a[0] + a[1]*exp(-a[2]*r) + a[3]*exp(-a[4]*r) + a[5]*exp(-a[6]*r))*inv_r;
}
