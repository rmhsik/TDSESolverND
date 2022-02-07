#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include "debug.h"
#include "wavefunction.h"

void WF::_geom_RZ(){
    _apply_mask = &WF::apply_mask_RZ;
    _pop = &WF::_pop_RZ;
    _pop_buf = &WF::_pop_buf_RZ;

    _acc_i = &WF::_acc_i_RZ;
    _acc_i_buf = &WF::_acc_i_buf_RZ;
    _acc_k = &WF::_acc_k_RZ;
    _acc_k_buf = &WF::_acc_k_buf_RZ;

    _dip_i = &WF::_dip_i_RZ;
    _dip_i_buf = &WF::_dip_i_buf_RZ;
    _dip_k = &WF::_dip_k_RZ;
    _dip_k_buf = &WF::_dip_k_buf_RZ;
}

void WF::apply_mask_RZ(cdouble *imask, cdouble *kmask){
    for(int i=0; i<_ni; i++){
        for(int k=0; k<_nk;k++)
            _wf[i*_nk + k] = _wf[i*_nk + k]*imask[i]*kmask[k];
    }
}


void WF::apply_mask_buf_RZ(cdouble *imask, cdouble *kmask,const int idx){
    for(int i=0; i<_ni; i++){
        for(int k=0; k<_nk;k++)
            _wf_buf[idx][i*_nk + k] *= imask[i]*kmask[k];
    }
}

cdouble WF::_acc_i_RZ(){
    cdouble sum = cdouble(0.0,0.0);
    //#pragma omp parallel for schedule(dynamics) reduction(+: sum)
    for(int i=0; i<_ni; i++){
        for(int k=0; k<_nk; k++){
            sum += 2*M_PI*_i[i]*conj(_wf[i*_nk + k])*(-1.0*_dV_i[i*_nk+k])*_wf[i*_nk+k]*_di*_dk;
        }
    }
    return sum;
}

cdouble WF::_acc_k_RZ(){
    cdouble sum = cdouble(0.0,0.0);
    //#pragma omp parallel for schedule(dynamics) reduction(+: sum)
    for(int i=0; i<_ni; i++){
        for(int k=0; k<_nk; k++){
            sum += 2*M_PI*_i[i]*conj(_wf[i*_nk + k])*(-1.0*_dV_k[i*_nk+k])*_wf[i*_nk+k]*_di*_dk;
        }
    }
    return sum;
}

void WF::_acc_i_buf_RZ(){
    #pragma omp parallel for schedule(dynamic)
    for(int n=0; n<_param->nt_diag;n++){
        cdouble sum = 0.0;
        for(int i=0;i<_ni;i++){
            for(int k=0;k<_nk;k++){
                sum += 2*M_PI*_i[i]*conj(_wf_buf[n][i*_nk + k])*(-1.0*_dV_i[i*_nk + k])*_wf_buf[n][i*_nk +k]*_di*_dk;
            }
        }
        _diag_buf[n] = sum;
    }
}

void WF::_acc_k_buf_RZ(){
    #pragma omp parallel for schedule(dynamic)
    for(int n=0; n<_param->nt_diag;n++){
        cdouble sum = 0.0;
        for(int i=0;i<_ni;i++){
            for(int k=0;k<_nk;k++){
                sum += 2*M_PI*_i[i]*conj(_wf_buf[n][i*_nk + k])*(-1.0*_dV_k[i*_nk + k])*_wf_buf[n][i*_nk +k]*_di*_dk;
            }
        }
        _diag_buf[n] = sum;
    }
} 


cdouble WF::_pop_RZ(double imin, double imax, double kmin, double kmax){
    cdouble sum = 0.0;
    int n_imin, n_imax;
    int n_kmin, n_kmax;
    n_imin = (imin-_i[0])/_di;
    n_imax = (imax-_i[0])/_di;
    n_kmin = (kmin-_k[0])/_dk;
    n_kmax = (kmax-_k[0])/_dk;
    for(int i=n_imin; i<n_imax;i++){
        for(int k=n_kmin;k<n_kmax;k++){
            sum += 2*M_PI*_i[i]*conj(_wf[i*_nk+k])*_wf[i*_nk+k]*_di*_dk;
        }
    }
    return sum;
}

void WF::_pop_buf_RZ(double imin, double imax, double kmin, double kmax){
    int n_imin, n_imax;
    int n_kmin, n_kmax;
    n_imin = (imin-_i[0])/_di;
    n_imax = (imax-_i[0])/_di;
    n_kmin = (kmin-_k[0])/_dk;
    n_kmax = (kmax-_k[0])/_dk;
    #pragma omp parallel for schedule(dynamic)
    for(int n=0;n<_param->nt_diag;n++){
        cdouble sum = 0.0;
        for(int i=n_imin; i<n_imax; i++){
            for(int k=n_kmin; k<n_kmax; k++){
                sum += 2*M_PI*_i[i]*conj(_wf_buf[n][i*_nk + k])*_wf_buf[n][i*_nk + k]*_di*_dk;
            }
        }
        _diag_buf[n] = sum;
    }
}

cdouble WF::_dip_i_RZ(){
    cdouble sum = 0.0;
    for(int i=0; i<_ni;i++){
        for(int k=0; k<_nk; k++){
            sum+=2*M_PI*_i[i]*conj(_wf[i*_nk+k])*_i[i]*_wf[i*_nk+k]*_di*_dk;
        }
    }
   return sum;
}

cdouble WF::_dip_k_RZ(){
    cdouble sum = 0.0;
    for(int i=0; i<_ni;i++){
        for(int k=0; k<_nk; k++){
            sum+=2*M_PI*_i[i]*conj(_wf[i*_nk+k])*_k[i]*_wf[i*_nk+k]*_di*_dk;
        }
    }
   return sum;
}

void WF::_dip_i_buf_RZ(){
    #pragma omp parallel for schedule(dynamic)
    for(int n=0; n<_param->nt_diag; n++){
        cdouble sum = 0.0;
        for(int i=0; i<_ni; i++){
            for(int k=0; k<_nk; k++){
                sum += 2*M_PI*_i[i]*conj(_wf_buf[n][i*_nk + k])*_i[k]*_wf_buf[n][i*_nk + k]*_di*_dk;
            }
        }
        _diag_buf[n] = sum;
    }
}

void WF::_dip_k_buf_RZ(){
    #pragma omp parallel for schedule(dynamic)
    for(int n=0; n<_param->nt_diag; n++){
        cdouble sum = 0.0;
        for(int i=0; i<_ni; i++){
            for(int k=0; k<_nk; k++){
                sum += 2*M_PI*_i[i]*conj(_wf_buf[n][i*_nk + k])*_k[k]*_wf_buf[n][i*_nk + k]*_di*_dk;
            }
        }
        _diag_buf[n] = sum;
    }

}

