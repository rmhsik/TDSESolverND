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
        Field *Afield_j;
        Field *Afield_k;
        Field *Bfield_i;
        Field *Bfield_j;
        Field *Bfield_k;

        double *_i, *_j, *_k, *_t;
        double _di, _dj, _dk, _dt;
        int _ni, _nj, _nk, _nt;
             
        cdouble *_Mk_du, *_Mk_d, *_Mk_dl;
        cdouble *_Mpk_du, *_Mpk_d, *_Mpk_dl;
        cdouble *_Mj_du, *_Mj_d, *_Mj_dl;
        cdouble *_Mpj_du, *_Mpj_d, *_Mpj_dl;
        cdouble *_Mi_du, *_Mi_d, *_Mi_dl;
        cdouble *_Mpi_du, *_Mpi_d, *_Mpi_dl;
        cdouble *_lhs_i, *_res_i;
        cdouble *_lhs_k, *_res_k; 
        cdouble *_lhs_j, *_res_j;
        void _allocate_X();
        void _allocate_XZ();
        void _allocate_RZ();
        void _allocate_XYZ();

        void tridot(cdouble* aa, cdouble *bb, cdouble* cc, cdouble* vec, cdouble* out, const int n);
        void tdma(cdouble* aa, cdouble *bb, cdouble* cc, cdouble* dd, cdouble* out,  const int n);

        cdouble (*_potential)(double i, double j, double k, double ti, Hamiltonian *ham);
        cdouble _potential_fn(double i, double j, double k, double ti);

    public:
        friend cdouble potential(double i, double j, double k, double ti, Hamiltonian *ham);
        friend cdouble potential_X(double i, double j, double k, double ti, Hamiltonian *ham);
        friend cdouble potential_XZ(double i, double j, double k, double ti, Hamiltonian *ham);
        friend cdouble potential_RZ(double i, double j, double k, double ti, Hamiltonian *ham);
        Hamiltonian();
        Hamiltonian(Parameters *param);
        void set_geometry(double *i, double *j, double *k, double *t, const double di, const double dj, const double dk, const double dt);
        void set_fields(Field* field1, Field* field2, Field* field3, Field* field4, Field* field5, Field* field6);
        cdouble dpotential_i(double i, double j, double k);
        cdouble dpotential_k(double i, double j, double k); 
        cdouble dpotential_j(double i, double j, double k);

        void (Hamiltonian::*step_i)(cdouble*, const int, const int, const int, const int, const int);
        void (Hamiltonian::*step_j)(cdouble*, const int, const int, const int, const int, const int);
        void (Hamiltonian::*step_k)(cdouble*, const int, const int, const int, const int, const int);
        cdouble (Hamiltonian::*ener)(cdouble ***psi);

        // Hamiltonian for X 
        void step_i_X(cdouble *psi, const int j, const int k, const int ti, const int imag, const int id_thread);
        void step_j_X(cdouble *psi, const int i, const int k, const int ti, const int imag, const int id_thread);
        void step_k_X(cdouble *psi, const int i, const int j, const int ti, const int imag, const int id_thread);
        cdouble ener_X(cdouble*** psi);

        // Hamiltonian for XZ 
        void step_i_XZ(cdouble *psi, const int j, const int k, const int ti, const int imag, const int id_thread);
        void step_j_XZ(cdouble *psi, const int i, const int k, const int ti, const int imag, const int id_thread);
        void step_k_XZ(cdouble *psi, const int i, const int j, const int ti, const int imag, const int id_thread);
        cdouble ener_XZ(cdouble*** psi);

        // Hamiltonian for RZ 
        void step_i_RZ(cdouble *psi, const int j, const int k, const int ti, const int imag, const int id_thread);
        void step_j_RZ(cdouble *psi, const int i, const int k, const int ti, const int imag, const int id_thread);
        void step_k_RZ(cdouble *psi, const int i, const int j, const int ti, const int imag, const int id_thread);
        cdouble ener_RZ(cdouble*** psi);

        // Hamiltonian for XYZ
        void step_i_XYZ(cdouble *psi, const int j, const int k, const int ti, const int imag, const int id_thread);
        void step_j_XYZ(cdouble *psi, const int i, const int k, const int ti, const int imag, const int id_thread);
        void step_k_XYZ(cdouble *psi, const int i, const int j, const int ti, const int imag, const int id_thread);
        cdouble ener_XYZ(cdouble*** psi);

        ~Hamiltonian();        
};

#endif
