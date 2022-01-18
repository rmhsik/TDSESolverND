#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H


#define cdouble std::complex<double>
#define I std::complex<double>(0.0,1.0)
#include <complex>
#include <string>
#include "parameters.h"

class WF{
    private:
        cdouble *_wf;
        cdouble *_wf_buf;
        cdouble *_diag_buf;
        cdouble *_row;
        cdouble *_col;
        cdouble *_dV;
        Parameters _param;
        int _ni, _nk;
        double *_i, *_k, _di, _dk;
        void (WF::*_apply_mask)(cdouble*,cdouble*);
    public:
        WF();
        WF(Parameters param);
        void set_geometry(double *i, double *k, const double di, const double dk);
        void gaussian(double i0, double k0, double sigma);
        void exponential(double i0, double k0, double sigma);
        cdouble* get();
        cdouble* row(int i);
        cdouble* col(int i); 
        void set(cdouble* arr);
        void set_row(cdouble* row, int k);
        void set_col(cdouble* col, int i);
        void set_to_buf(const int idx);
        void set_col_buf(cdouble* col, const int i, const int idx);
        void set_row_buf(cdouble* row, const int k, const int idx);
        void get_from_buf(cdouble* arr, const int idx);
        cdouble* get_buf();
        cdouble* get_diag_buf();
        cdouble norm();
        void set_dpotential(cdouble *dV);
        void apply_mask(cdouble* imask, cdouble *kmask);
        void apply_mask_X(cdouble* imask, cdouble *kmask);
        void apply_mask_RZ(cdouble* imask, cdouble *kmask);
        void apply_mask_buf_RZ(cdouble* imask, cdouble* kmask, const int idx);
        void save_wf(std::string path);
        void save_wf2(std::string path);
        cdouble operator()(int i, int j);
        void operator/= (cdouble val);        

        cdouble dipole();
        cdouble acc();
        
        void dipole_buf();
        void acc_buf();

};



#endif
