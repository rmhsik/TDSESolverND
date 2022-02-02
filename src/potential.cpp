#include <cmath>
#include "potential.h"


cdouble potential( double i, double k, int ti, Hamiltonian *ham){
    return -1.0/sqrt(i*i + k*k);
}

cdouble potential_X( double i, double k, int ti, Hamiltonian *ham){
    return -1.0/sqrt(i*i + k*k + 2.0);
}
cdouble potential_XZ( double i, double k, int ti, Hamiltonian *ham){
    return -1.0/sqrt(i*i + k*k + 0.65);
}
cdouble potential_RZ( double i, double k, int ti, Hamiltonian *ham){
    return -1.0/sqrt(i*i + k*k );
}
