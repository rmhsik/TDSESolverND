#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H


#define cdouble std::complex<double>
#define I std::complex<double>(0.0,1.0)
#include <complex>
#include <string>
#include "parameters.h"

class WF{
    private:
        cdouble **_wf;
        cdouble ***_wf_buf;
        cdouble *_diag_buf;
        cdouble *_i_row;
        cdouble *_k_row;
        cdouble *_dV_i;
        cdouble *_dV_k;
        Parameters *_param;
        int _ni, _nk;
        double *_i, *_k, _di, _dk;
        void (WF::*_apply_mask)(cdouble*,cdouble*);

        void _geom_XZ();

    public:
        WF();
        WF(Parameters *param);
        void set_geometry(double *i, double *k, const double di, const double dk);
        void set_diagnostics();
        void gaussian(double i0, double k0, double sigma);
        void gaussian_anti(double i0, double k0, double sigma);
        void exponential(double i0, double k0, double sigma);
        cdouble** get();
        cdouble* i_row(int k);
        cdouble* k_row(int i);
        void set(cdouble** arr);
        void set_i_row(cdouble* i_row, int k);
        void set_i_row_mask(cdouble* i_row, cdouble *imask, int k);
        void set_k_row(cdouble* k_row, int i);
        void set_k_row_mask(cdouble* k_row, cdouble *kmask, int i);
        void get_i_row(cdouble* i_row, int k);
        void get_k_row(cdouble* k_row, int i);

        void set_to_buf(const int idx);
        void set_i_row_buf(cdouble* i_row, const int k, const int idx);
        void set_k_row_buf(cdouble* k_row, const int i, const int idx);
        void get_i_row_buf(cdouble* i_row, const int k, const int idx);
        void get_k_row_buf(cdouble* k_row, const int i, const int idx);


        void set_i_row_buf_mask(cdouble* i_row, cdouble* imask, const int k, const int idx);
        void set_k_row_buf_mask(cdouble* k_row, cdouble* kmask, const int i, const int idx);
        void get_from_buf(cdouble** arr, const int idx);

	void anti_sym_k();

        cdouble*** get_buf();

        cdouble norm();
	cdouble norm_buf(const int idx);
        void apply_mask(cdouble* imask, cdouble *kmask);
        void apply_mask_XZ(cdouble* imask, cdouble *kmask);
        void apply_mask_buf_XZ(cdouble* imask, cdouble* kmask, const int idx);
        cdouble operator()(int i, int k);
        void operator/= (cdouble val);        
	void save_wf2(std::string name);
	~WF();
};



#endif
