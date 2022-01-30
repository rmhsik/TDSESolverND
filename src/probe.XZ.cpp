#include "probe.XZ.h"

ProbeXZ::ProbeXZ(): Probe(){}

ProbeXZ::ProbeXZ(std::string def): Probe::Probe(def){}

void ProbeXZ::_acc_i(const int idx){
    cdouble* wf_buf = _wf->get_buf();
    cdouble* dV_i = _ham->get_dpotential_i();
    #pragma omp parallel for schedule(dynamic)
    for(int n=0;n<_nt_diag;n++){
        cdouble sum = 0.0; 
        for(int i=0; i<_ni; i++){
            for(int k=0; k< _nk; k++){
                sum += conj(wf_buf[n*_ni*_nk + i*_nk + k])*(-1.0*dV_i[i*_nk + k]) * wf_buf[n*_ni*_nk + i*_nk + k]*_di*_dk;
            }
        }
        _data[idx + n] = sum*_tempmask[idx+n];
    }
}

void ProbeXZ::_acc_k(const int idx){
    cdouble* wf_buf = _wf->get_buf();
    cdouble* dV_k = _ham->get_dpotential_k();
    #pragma omp parallel for schedule(dynamic)
    for(int n=0;n<_nt_diag;n++){
        cdouble sum = 0.0; 
        for(int i=0; i<_ni; i++){
            for(int k=0; k< _nk; k++){
                sum += conj(wf_buf[n*_ni*_nk + i*_nk + k])*(-1.0*dV_k[i*_nk + k]) * wf_buf[n*_ni*_nk + i*_nk + k]*_di*_dk;
            }
        }
        _data[idx + n] = sum*_tempmask[idx+n];
    }
}

void ProbeXZ::_dip_i(const int idx){
    int n_imin = (_int_imin - _i[0])/_di;
    int n_imax = (_int_imax - _i[0])/_di;
    int n_kmin = (_int_kmin - _k[0])/_dk;
    int n_kmax = (_int_kmax - _k[0])/_dk;

    cdouble* wf_buf = _wf->get_buf();
    #pragma omp parallel for schedule(dynamic)
    for(int n=0;n<_nt_diag;n++){
        cdouble sum = 0.0; 
        for(int i=n_imin; i<n_imax; i++){
            for(int k=n_kmin; k<n_kmax; k++){
                sum += conj(wf_buf[n*_ni*_nk + i*_nk + k])*(-1.0*_i[i]) * wf_buf[n*_ni*_nk + i*_nk + k]*_di*_dk;
            }
        }
        _data[idx + n] = sum*_tempmask[idx+n];
    }
}

void ProbeXZ::_dip_k(const int idx){
    int n_imin = (_int_imin - _i[0])/_di;
    int n_imax = (_int_imax - _i[0])/_di;
    int n_kmin = (_int_kmin - _k[0])/_dk;
    int n_kmax = (_int_kmax - _k[0])/_dk;

    cdouble* wf_buf = _wf->get_buf();
    #pragma omp parallel for schedule(dynamic)
    for(int n=0;n<_nt_diag;n++){
        cdouble sum = 0.0; 
        for(int i=n_imin; i<n_imax; i++){
            for(int k=n_kmin; k<n_kmax; k++){
                sum += conj(wf_buf[n*_ni*_nk + i*_nk + k])*(-1.0*_k[k]) * wf_buf[n*_ni*_nk + i*_nk + k]*_di*_dk;
            }
        }
        _data[idx + n] = sum*_tempmask[idx+n];
    }
}

void ProbeXZ::_pop(const int idx){
    int n_imin = (_int_imin - _i[0])/_di;
    int n_imax = (_int_imax - _i[0])/_di;
    int n_kmin = (_int_kmin - _k[0])/_dk;
    int n_kmax = (_int_kmax - _k[0])/_dk;

    cdouble* wf_buf = _wf->get_buf();
    #pragma omp parallel for schedule(dynamic)
    for(int n=0;n<_nt_diag;n++){
        cdouble sum = 0.0; 
        for(int i=n_imin; i<n_imax; i++){
            for(int k=n_kmin; k<n_kmax; k++){
                sum += conj(wf_buf[n*_ni*_nk + i*_nk + k])* wf_buf[n*_ni*_nk + i*_nk + k]*_di*_dk;
            }
        }
        _data[idx + n] = sum;
    }
}
