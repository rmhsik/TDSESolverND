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
        cdouble **_wf_buf;
        cdouble *_diag_buf;
        cdouble *_row;
        cdouble *_col;
        cdouble *_dV_i;
        cdouble *_dV_k;
        Parameters *_param;
        int _ni, _nk;
        double *_i, *_k, _di, _dk;
        void (WF::*_apply_mask)(cdouble*,cdouble*);

        cdouble (WF::*_pop)(double,double,double,double);
        void (WF::*_pop_buf)(double,double,double,double);
        cdouble (WF::*_acc_i)();
        cdouble (WF::*_acc_k)();
        cdouble (WF::*_dip_i)();
        cdouble (WF::*_dip_k)();
        void (WF::*_acc_i_buf)();
        void (WF::*_acc_k_buf)();
        void (WF::*_dip_i_buf)();
        void (WF::*_dip_k_buf)();
        
	void _geom_X();
	void _geom_XZ();
	void _geom_RZ();

        cdouble _pop_X(double imin, double imax, double kmin, double kmax);
        cdouble _pop_XZ(double imin, double imax, double kmin, double kmax);
        cdouble _pop_RZ(double imin, double imax, double kmin, double kmax);
        cdouble _pop_0(double imin, double imax, double kmin, double kmax);
        cdouble _acc_i_X();
        cdouble _acc_k_X();
        cdouble _acc_i_XZ();
        cdouble _acc_k_XZ();
        cdouble _acc_i_RZ();
        cdouble _acc_k_RZ();
        cdouble _acc_0();
        cdouble _dip_i_X();
        cdouble _dip_k_X();
        cdouble _dip_i_XZ();
        cdouble _dip_k_XZ();
        cdouble _dip_i_RZ();
        cdouble _dip_k_RZ();
        cdouble _dip_0();

        void _acc_i_buf_X();
        void _acc_k_buf_X();
        void _acc_i_buf_XZ();
        void _acc_k_buf_XZ();
        void _acc_i_buf_RZ();
        void _acc_k_buf_RZ();
        void _acc_buf_0();
        void _pop_buf_X(double imin, double imax, double kmin, double kmax);
        void _pop_buf_XZ(double imin, double imax, double kmin, double kmax);
        void _pop_buf_RZ(double imin, double imax, double kmin, double kmax);
        void _pop_buf_0(double imin, double imax, double kmin, double kmax);
        void _dip_i_buf_X();
        void _dip_k_buf_X();
        void _dip_i_buf_XZ();
        void _dip_k_buf_XZ();
        void _dip_i_buf_RZ();
        void _dip_k_buf_RZ();
        void _dip_buf_0();

    public:
        WF();
        WF(Parameters *param);
        void set_geometry(double *i, double *k, const double di, const double dk);
        void set_diagnostics();
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
        void set_col_buf_mask(cdouble* col, cdouble* kmask, const int i, const int idx);
        void set_row_buf_mask(cdouble* row, cdouble* imask, const int k, const int idx);
        void get_from_buf(cdouble* arr, const int idx);
        cdouble** get_buf();
        cdouble* get_diag_buf();
        cdouble norm();
        void apply_mask(cdouble* imask, cdouble *kmask);
        void apply_mask_X(cdouble* imask, cdouble *kmask);
        void apply_mask_RZ(cdouble* imask, cdouble *kmask);
        void apply_mask_buf_RZ(cdouble* imask, cdouble* kmask, const int idx);
        void apply_mask_XZ(cdouble* imask, cdouble *kmask);
        void apply_mask_buf_XZ(cdouble* imask, cdouble* kmask, const int idx);
        void save_wf(std::string path);
        void save_wf2(std::string path);
        cdouble operator()(int i, int j);
        void operator/= (cdouble val);        

        cdouble dip_i();
        cdouble dip_k();
        cdouble acc_i();
        cdouble acc_k();
        cdouble pop(double imin, double imax, double kmin, double kmax);
        
        void dip_i_buf();
        void dip_k_buf();
        void acc_i_buf();
        void acc_k_buf();
        void pop_buf(double imin, double imax, double kmin, double kmax);

	~WF();
};



#endif
