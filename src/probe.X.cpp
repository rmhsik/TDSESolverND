#include "probe.X.h"

ProbeX::ProbeX(): Probe::Probe(){}

ProbeX::ProbeX(std::string def): Probe::Probe(def){}

void ProbeX::_acc_i(const int idx){
    cdouble* wf_buf = _wf->get_buf();
    //cdouble* dV_i = _ham->get_dpotential_i();
    for(int n=0; n<_param->nt_diag;n++){
        cdouble sum = 0.0;
        cdouble dV_i=0.0;
        for(int i=0; i<_ni; i++){
            dV_i = _ham->dpotential_i(_i[i],0);
            sum += conj(wf_buf[n*_ni*_nk + i*_nk + 0])*(-1.0*dV_i)*wf_buf[n*_ni*_nk + i*_nk +0]*_di;
        }
        _data[idx + n] = sum;
    }
}

void ProbeX::_acc_k(const int idx){}

void ProbeX::_dip_i(const int idx){
    int n_imin = (_int_imin - _i[0])/_di;
    int n_imax = (_int_imax - _i[0])/_di;
    int n_kmin = (_int_kmin - _k[0])/_dk;
    int n_kmax = (_int_kmax - _k[0])/_dk;


    cdouble* wf_buf = _wf->get_buf();
    for(int n=0; n<_param->nt_diag;n++){
        cdouble sum = 0.0;
        for(int i=n_imin; i<n_imax; i++){
            sum += conj(wf_buf[n*_ni*_nk + i*_nk + 0])*(-1.0*_i[i])*wf_buf[n*_ni*_nk + i*_nk +0]*_di;
        }
        _data[idx + n] = sum;
    }
}

void ProbeX::_dip_k(const int idx){}

void ProbeX::_pop(const int idx){
    int n_imin = (_int_imin - _i[0])/_di;
    int n_imax = (_int_imax - _i[0])/_di;
    int n_kmin = (_int_kmin - _k[0])/_dk;
    int n_kmax = (_int_kmax - _k[0])/_dk;


    cdouble* wf_buf = _wf->get_buf();
    for(int n=0; n<_param->nt_diag;n++){
        cdouble sum = 0.0;
        for(int i=n_imin; i<n_imax; i++){
            sum += conj(wf_buf[n*_ni*_nk + i*_nk + 0])*wf_buf[n*_ni*_nk + i*_nk +0]*_di;
        }
        _data[idx+n] = sum;
    }
}
