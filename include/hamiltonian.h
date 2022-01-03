#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "wavefunction.h"
class Hamiltonian{
    private:
        Parameters _param;
        WF *_wf;
        cdouble **_potential;
    public:
        Hamiltonian();
        Hamiltonian(Parameters param, WF *wf, double *i, double *k, const double di, const double dk);
        

};

#endif
