#include <cmath>
#include "potential.h"

#define A 2.0

cdouble potential( double i, double k, int ti, Hamiltonian *ham){
    return -1.0/sqrt(i*i + k*k + A);
}
