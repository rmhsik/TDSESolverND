#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#define cdouble std::complex<double>
#define I std::complex<double>(0.0,1.0)
#include <complex>
#include <string>
#include "parameters.h"

class Hamiltonian{
    private:
        Parameters _param;
        cdouble **_potential;
        cdouble **_dpotential;
        double *_i, *_k;
        double _di, _dk;
        int _ni, _nk;
        
        void potential();
        void dpotential();
        void tridot(cdouble* aa, cdouble *bb, cdouble* cc, cdouble* vec, cdouble* out, const int n);
        void tdma(cdouble* aa, cdouble *bb, cdouble* cc, cdouble* dd, cdouble* out,  const int n);
    public:
        Hamiltonian();
        Hamiltonian(Parameters param);
        void set_geometry(double *i, double *k, const double di, const double dk);
        void set_potential();
        void set_dpotential();
        cdouble **get_potential();
        cdouble **get_dpotential();
        void step_i(cdouble *psi, double afield_i, double bfield_i, const int imag);
        void step_k(cdouble *psi);
        cdouble ener(cdouble* psi);
};

#endif
