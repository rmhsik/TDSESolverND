#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H


#define cdouble std::complex<double>
#define I std::complex<double>(0.0,1.0)
#include <complex>
#include <string>
#include "parameters.h"

class WF{
    private:
        cdouble ***_wf;
	cdouble ***_wf_0;
        cdouble ****_wf_buf;
        cdouble *_diag_buf;
        cdouble *_i_row;
        cdouble *_j_row;
        cdouble *_k_row;
        cdouble *_dV_i;
        cdouble *_dV_k;
        Parameters *_param;
        int _ni, _nj ,_nk;
        double *_i, *_j, *_k, _di, _dj, _dk;
        void (WF::*_apply_mask)(cdouble*,cdouble*,cdouble*);

        void _geom_X();
        void _geom_XZ();
        void _geom_RZ();
        void _geom_XYZ();

    public:
        WF();
        WF(Parameters *param);
        void set_geometry(double *i, double *j, double *k, const double di, const double dj, const double dk);
        void set_diagnostics();
        void gaussian(double i0, double j0, double k0, double sigma);
        void exponential(double i0, double j0, double k0, double sigma);
        cdouble*** get();
        cdouble* i_row(int j, int k);
        cdouble* j_row(int i, int k); 
        cdouble* k_row(int i, int j);
        void set(cdouble*** arr);
        void set_i_row(cdouble* i_row, int j, int k);
        void set_i_row_mask(cdouble* i_row, cdouble *imask, int j, int k);
        void set_j_row(cdouble* j_row, int i, int k);
        void set_j_row_mask(cdouble* j_row, cdouble *jmask, int i, int k);
        void set_k_row(cdouble* k_row, int i, int j);
        void set_k_row_mask(cdouble* k_row, cdouble *kmask, int i, int j);
        void get_i_row(cdouble* i_row, int j, int k);
        void get_j_row(cdouble* j_row, int i, int k);
        void get_k_row(cdouble* k_row, int i, int j);

        void set_to_buf(const int idx);
        void set_i_row_buf(cdouble* i_row, const int j, const int k, const int idx);
        void set_j_row_buf(cdouble* j_row, const int i, const int k, const int idx);
        void set_k_row_buf(cdouble* k_row, const int i, const int j, const int idx);
        void get_i_row_buf(cdouble* i_row, const int j, const int k, const int idx);
        void get_j_row_buf(cdouble* j_row, const int i, const int k, const int idx);
        void get_k_row_buf(cdouble* k_row, const int i, const int j, const int idx);


        void set_i_row_buf_mask(cdouble* i_row, cdouble* imask, const int j, const int k, const int idx);
        void set_j_row_buf_mask(cdouble* j_row, cdouble* jmask, const int i, const int k, const int idx);
        void set_k_row_buf_mask(cdouble* k_row, cdouble* kmask, const int i, const int j, const int idx);
        void get_from_buf(cdouble*** arr, const int idx);

        cdouble**** get_buf();
        cdouble* get_diag_buf();

	void set_to_ground();
	cdouble ***get_ground();

        cdouble norm();
        void apply_mask(cdouble* imask, cdouble *jmask, cdouble *kmask);
        void apply_mask_X(cdouble* imask, cdouble *jmask, cdouble *kmask);
        void apply_mask_RZ(cdouble* imask, cdouble *jmask, cdouble *kmask);
        void apply_mask_buf_RZ(cdouble* imask, cdouble *jmask, cdouble* kmask, const int idx);
        void apply_mask_XZ(cdouble* imask, cdouble *jmask, cdouble *kmask);
        void apply_mask_buf_XZ(cdouble* imask, cdouble *jmask, cdouble* kmask, const int idx);
        void apply_mask_XYZ(cdouble* imask, cdouble *jmask, cdouble *kmask);
        void apply_mask_buf_XYZ(cdouble* imask, cdouble *jmask, cdouble* kmask, const int idx);
        cdouble operator()(int i, int j, int k);
        void operator/= (cdouble val);        
	~WF();
};



#endif
