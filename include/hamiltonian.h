#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#define cdouble std::complex<double>
#define I std::complex<double>(0.0,1.0)
#include <complex>
#include <string>
#include "parameters.h"
#include "fields.h"

class Hamiltonian{
    private:
        Parameters *_param;

        Field *Afield_i;
        Field *Afield_k;

        double *_i, *_k, *_t;
        double _di, _dk, _dt;
        int _ni, _nk, _nt;
             
        cdouble *_Mk_du, *_Mk_d, *_Mk_dl;
        cdouble *_Mpk_du, *_Mpk_d, *_Mpk_dl;
        cdouble *_Mi_du, *_Mi_d, *_Mi_dl;
        cdouble *_Mpi_du, *_Mpi_d, *_Mpi_dl;
        cdouble *_lhs_i, *_res_i;
        cdouble *_lhs_k, *_res_k; 
        void _allocate_XZ();

        void tridot(cdouble* aa, cdouble *bb, cdouble* cc, cdouble* vec, cdouble* out, const int n);
        void tdma(cdouble* aa, cdouble *bb, cdouble* cc, cdouble* dd, cdouble* out,  const int n);

        cdouble (*_potential)(double i, double k, double ti, Hamiltonian *ham);
        cdouble _potential_fn(double i, double k, double ti);

    public:
        friend cdouble potential(double i, double k, double ti, Hamiltonian *ham);
        friend cdouble potential_XZ(double i, double k, double ti, Hamiltonian *ham);
        Hamiltonian();
        Hamiltonian(Parameters *param);
        void set_geometry(double *i, double *k, double *t, const double di, const double dk, const double dt);
        void set_fields(Field* field1, Field* field2);
        cdouble dpotential_i(double i, double k);
        cdouble dpotential_k(double i, double k); 

        void (Hamiltonian::*step_i)(cdouble*, const int, const int, const int, const int);
        void (Hamiltonian::*step_k)(cdouble*, const int, const int, const int, const int);
        cdouble (Hamiltonian::*ener)(cdouble **psi);

        // Hamiltonian for XZ 
        void step_i_XZ(cdouble *psi, const int k, const int ti, const int imag, const int id_thread);
        void step_k_XZ(cdouble *psi, const int i, const int ti, const int imag, const int id_thread);
        cdouble ener_XZ(cdouble** psi);

        ~Hamiltonian();        
};

#endif
