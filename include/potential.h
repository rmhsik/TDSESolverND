#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <cmath>
#include "hamiltonian.h"

cdouble potential(double i, double k, int ti, Hamiltonian *ham);
cdouble potential_X(double i, double k, int ti, Hamiltonian *ham);
cdouble potential_XZ(double i, double k, int ti, Hamiltonian *ham);
cdouble potential_RZ(double i, double k, int ti, Hamiltonian *ham);
#endif

