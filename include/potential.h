#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <cmath>
#include "hamiltonian.h"

cdouble potential(double i, double j, double k, double ti, Hamiltonian *ham);
cdouble potential_X(double i, double j, double k, double t, Hamiltonian *ham);
cdouble potential_XZ(double i, double j, double k, double t, Hamiltonian *ham);
cdouble potential_RZ(double i, double j, double k, double t, Hamiltonian *ham);
#endif

