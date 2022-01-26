#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include "debug.h"
#include "wavefunction.h"

void WF::_geom_X(){
    _apply_mask = &WF::apply_mask_X;
    _pop = &WF::_pop_X;
    _pop_buf = &WF::_pop_buf_X;

    _acc_i = &WF::_acc_i_X;
    _acc_i_buf = &WF::_acc_i_buf_X;
    _acc_k = &WF::_acc_k_X;
    _acc_k_buf = &WF::_acc_k_buf_X;
    
    _dip_i = &WF::_dip_i_X;
    _dip_i_buf = &WF::_dip_i_buf_X;
    _dip_k = &WF::_dip_k_X;
    _dip_k_buf = &WF::_dip_k_buf_X;

}

void WF::apply_mask_X(cdouble *imask, cdouble *kmask){
    for(int i=0; i<_ni; i++){
        _wf[i*_nk + 0] = _wf[i*_nk + 0]*imask[i];
    }
}

cdouble WF::_acc_i_X(){return 0.0;}

cdouble WF::_acc_k_X(){
    cdouble sum = cdouble(0.0,0.0);
    for(int i=0; i<_ni; i++){
        sum+= conj(_wf[i*_nk + 0])*(-1.0*_dV_k[i*_nk + 0])*_wf[i*_nk + 0]*_di;
    }
    return sum;
}

void WF::_acc_i_buf_X(){}

void WF::_acc_k_buf_X(){
    for(int n=0; n<_param->nt_diag;n++){
        cdouble sum = 0.0;
        for(int i=0; i<_ni; i++){
            sum += conj(_wf_buf[n*_ni*_nk + i*_nk + 0])*(-1.0*_dV_k[i*_nk + 0])*_wf_buf[n*_ni*_nk + i*_nk +0]*_di;
        }
        _diag_buf[n] = sum;
    }
}

cdouble WF::_pop_X(double imin, double imax, double kmin, double kmax){
    cdouble sum = 0.0;
    int n_imin, n_imax;
    n_imin = (imin-_i[0])/_di;
    n_imax = (imax-_i[0])/_di;
    for(int i=n_imin; i<n_imax;i++){
        for(int k=0;k<1;k++){
            sum += conj(_wf[i*_nk+k])*_wf[i*_nk+k]*_di;
        }
    }
    return sum;
}

void WF::_pop_buf_X(double imin, double imax, double kmin, double kmax){
    int n_imin, n_imax;
    int n_kmin, n_kmax;
    n_imin = (imin-_i[0])/_di;
    n_imax = (imax-_i[0])/_di;
    n_kmin = (kmin-_k[0])/_dk;
    n_kmax = (kmax-_k[0])/_dk;
    for(int n=0;n<_param->nt_diag;n++){
        cdouble sum = 0.0;
        for(int i=n_imin; i<n_imax; i++){
            for(int k=0; k<1; k++){
                sum += conj(_wf_buf[n*_nk*_ni + i*_nk + k])*_wf_buf[n*_nk*_ni + i*_nk + k]*_di;
            }
        }
        _diag_buf[n] = sum;
    }
}

cdouble WF::_dip_i_X(){
    cdouble sum = 0.0;
    for(int i=0; i<_ni; i++){
        sum += conj(_wf[i*_nk +0])*_i[i]*_wf[i*_nk + 0]*_di;
    }
    return sum;
}

cdouble WF::_dip_k_X(){
    return 0.0;
}

void WF::_dip_i_buf_X(){
    for(int n=0; n<_param->nt_diag; n++){
        cdouble sum = 0.0;
        for(int i=0;i<_ni;i++){
            sum += conj(_wf_buf[n*_ni*_nk + i*_nk + 0])*_i[i]*_wf_buf[n*_ni*_nk + i*_nk + 0]*_di;
        }
        _diag_buf[n] = sum;
    }
}

void WF::_dip_k_buf_X(){}
