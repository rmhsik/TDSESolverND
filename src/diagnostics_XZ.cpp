#include <cmath>
#include "diagnostics.h"

DiagnosticsXZ::DiagnosticsXZ():Diagnostics(){}

DiagnosticsXZ::DiagnosticsXZ(Parameters *param):Diagnostics(param){
}

void DiagnosticsXZ::_acc_i_buf(const int idx){
    cdouble* wf_buf = _wf->get_buf();
    cdouble* dV_i = _ham->get_dpotential_i();

    for(int n=0;n<_param->nt_diag;n++){
        cdouble sum = 0.0; 
        for(int i=0; i<_ni; i++){
            for(int k=0; k< _nk; k++){
                sum += conj(wf_buf[n*_ni*_nk + i*_nk + k])*(-1.0*dV_i[i*_nk + k]) * wf_buf[n*_ni*_nk + i*_nk + k]*_di*_dk;
            }
        }
        _acc_i[idx + n] = sum;
    }
}

void DiagnosticsXZ::_acc_k_buf(const int idx){
    cdouble* wf_buf = _wf->get_buf();
    cdouble* dV_k = _ham->get_dpotential_k();

    for(int n=0;n<_param->nt_diag;n++){
        cdouble sum = 0.0; 
        for(int i=0; i<_ni; i++){
            for(int k=0; k< _nk; k++){
                sum += conj(wf_buf[n*_ni*_nk + i*_nk + k])*(-1.0*dV_k[i*_nk + k]) * wf_buf[n*_ni*_nk + i*_nk + k]*_di*_dk;
            }
        }
        _acc_k[idx + n] = sum;
    }
}
void DiagnosticsXZ::_dip_i_buf(const int idx){
    cdouble* wf_buf = _wf->get_buf();

    for(int n=0;n<_param->nt_diag;n++){
        cdouble sum = 0.0; 
        for(int i=0; i<_ni; i++){
            for(int k=0; k< _nk; k++){
                sum += conj(wf_buf[n*_ni*_nk + i*_nk + k])*_i[i] * wf_buf[n*_ni*_nk + i*_nk + k]*_di*_dk;
            }
        }
        _dip_i[idx + n] = sum;
    }
}

void DiagnosticsXZ::_dip_k_buf(const int idx){
    cdouble* wf_buf = _wf->get_buf();

    for(int n=0;n<_param->nt_diag;n++){
        cdouble sum = 0.0; 
        for(int i=0; i<_ni; i++){
            for(int k=0; k< _nk; k++){
                sum += conj(wf_buf[n*_ni*_nk + i*_nk + k])*_k[i] * wf_buf[n*_ni*_nk + i*_nk + k]*_di*_dk;
            }
        }
        _dip_k[idx + n] = sum;
    }
}

void DiagnosticsXZ::_pop_buf(const int idx){
    cdouble* wf_buf = _wf->get_buf();
    int n_imin, n_imax;
    int n_kmin, n_kmax;
    n_imin = (_param->pop_imin - _i[0])/_di;
    n_imax = (_param->pop_imax - _i[0])/_di;
    n_kmin = (_param->kmin - _k[0])/_dk;
    n_kmax = (_param->kmax - _k[0])/_dk;

    for(int n=0; n<_param->nt_diag; n++){
        cdouble sum = 0.0;
        for(int i=n_imin; i<n_imax; i++){
            for(int k=n_kmin; k<n_kmax; k++){
                sum += conj(wf_buf[n*_ni*_nk + i*_nk + k])*wf_buf[n*_ni*_nk + i*_nk + k];
            }
        }
        _pop[idx + n] = sum;
    }
}
