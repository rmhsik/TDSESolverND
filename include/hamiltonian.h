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

        void (Hamiltonian::*step_i)(cdouble*, double, double, const int, const int);
        void (Hamiltonian::*step_k)(cdouble*, double, double, const int, const int);
        cdouble (Hamiltonian::*ener)(cdouble **psi);

        // Hamiltonian for X 
        void step_i_X(cdouble *psi_row, double afield_i, double bfield_i, const int j ,const int imag);
        void step_k_X(cdouble *psi_col, double afield_k, double bfield_k, const int i, const int imag);
        cdouble ener_X(cdouble** psi);

        // Hamiltonian for RZ 
        void step_i_RZ(cdouble *psi_row, double afield_i, double bfield_i, const int j ,const int imag);
        void step_k_RZ(cdouble *psi_col, double afield_k, double bfield_k, const int i, const int imag);
        cdouble ener_RZ(cdouble** psi);
};

#endif
