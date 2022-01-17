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
        cdouble *_potential;
        cdouble *_dpotential;
        double *_i, *_k;
        double _di, _dk;
        int _ni, _nk;
        
        cdouble *_Mk_du, *_Mk_d, *_Mk_dl;
        cdouble *_Mpk_du, *_Mpk_d, *_Mpk_dl;
        cdouble *_Mi_du, *_Mi_d, *_Mi_dl;
        cdouble *_Mpi_du, *_Mpi_d, *_Mpi_dl;
        cdouble *_lhs_i, *_res_i;
        cdouble *_lhs_k, *_res_k; 

        void potential_X();
        void dpotential_X();
        void potential_RZ();
        void dpotential_RZ();

        void tridot(cdouble* aa, cdouble *bb, cdouble* cc, cdouble* vec, cdouble* out, const int n);
        void tdma(cdouble* aa, cdouble *bb, cdouble* cc, cdouble* dd, cdouble* out,  const int n);
    public:
        Hamiltonian();
        Hamiltonian(Parameters param);
        void set_geometry(double *i, double *k, const double di, const double dk);
        void set_potential();
        void set_dpotential();
        cdouble *get_potential();
        cdouble *get_dpotential();

        void (Hamiltonian::*step_i)(cdouble*, double, double, const int, const int, const int);
        void (Hamiltonian::*step_k)(cdouble*, double, double, const int, const int, const int);
        cdouble (Hamiltonian::*ener)(cdouble *psi);

        // Hamiltonian for X 
        void step_i_X(cdouble *psi, double afield_i, double bfield_i, const int j ,const int imag, const int id_thread);
        void step_k_X(cdouble *psi, double afield_k, double bfield_k, const int i, const int imag, const int id_thread);
        cdouble ener_X(cdouble* psi);

        // Hamiltonian for RZ 
        void step_i_RZ(cdouble *psi, double afield_i, double bfield_k, const int j ,const int imag, const int id_thread);
        void step_k_RZ(cdouble *psi, double afield_k, double bfield_k, const int i, const int imag, const int id_thread);
        cdouble ener_RZ(cdouble* psi);
};

#endif
