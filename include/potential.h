#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <cmath>
#include "hamiltonian.h"

cdouble potential(double i, double k, double ti, Hamiltonian *ham);
cdouble potential_hydrogen_XZ(double i, double k, double t, Hamiltonian *ham);
cdouble potential_helium_XZ(double i, double k, double t, Hamiltonian *ham);
#endif

