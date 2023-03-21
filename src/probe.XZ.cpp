#include "probe.XZ.h"
#include <iostream>
#include <omp.h>

ProbeXZ::ProbeXZ(): Probe(){}

ProbeXZ::ProbeXZ(std::string def): Probe::Probe(def){}

void ProbeXZ::_acc_i(const int idx){
    cdouble** wf = _wf->get();
    cdouble *sum_tmp;
    sum_tmp = new cdouble[_param->n_threads];
    cdouble **dV_i;
    dV_i = _ham->dpotential_i();
    for(int i=0; i<_param->n_threads; i++){
        sum_tmp[i] = cdouble(0.0,0.0);
    }
    
	#pragma omp parallel for schedule(static) collapse(2)
    for(int i=0; i<_ni; i++){
        for(int k=0; k< _nk; k++){
            int idthread = omp_get_thread_num();
            //cdouble dV_i = _ham->dpotential_i(i,k);
            sum_tmp[idthread] += conj(wf[i][k])*(-1.0*dV_i[i][k]) * wf[i][k]*_di*_dk;
        }
    }
    
    cdouble sum = cdouble(0.0,0.0);
    for(int i=0; i<_param->n_threads; i++){
        sum += sum_tmp[i];
    }

    _data[idx] = sum*_tempmask[idx];
    std::cout<<"acc_i: "<<_data[idx]<<std::endl;
}

void ProbeXZ::_acc_k(const int idx){
    cdouble** wf = _wf->get();
    cdouble *sum_tmp;
    sum_tmp = new cdouble[_param->n_threads];
    cdouble **dV_k;
    dV_k  = _ham->dpotential_k();
    for(int i=0; i<_param->n_threads; i++){
        sum_tmp[i] = cdouble(0.0,0.0);
    }

	#pragma omp parallel for schedule(static) collapse(2)
    for(int i=0; i<_ni; i++){
        for(int k=0; k< _nk; k++){
            int idthread = omp_get_thread_num();
            //cdouble dV_k = _ham->dpotential_k(i,k);
            sum_tmp[idthread] += conj(wf[i][k])*(-1.0*dV_k[i][k]) * wf[i][k]*_di*_dk;
        }
    }

    cdouble sum = cdouble(0.0,0.0);
    for(int i=0; i<_param->n_threads; i++){
        sum += sum_tmp[i];
    }

    _data[idx] = sum*_tempmask[idx];
    std::cout<<"acc_k: "<<_data[idx]<<std::endl;
    
}

void ProbeXZ::_dip_i(const int idx){
    int n_imin = (_int_imin - _i[0])/_di;
    int n_imax = (_int_imax - _i[0])/_di;
    int n_kmin = (_int_kmin - _k[0])/_dk;
    int n_kmax = (_int_kmax - _k[0])/_dk;

    cdouble*** wf_buf = _wf->get_buf();
    for(int n=0;n<_nt_diag;n++){
        cdouble sum = 0.0; 
	#pragma omp parallel for schedule(dynamic) collapse(1)
        for(int i=n_imin; i<n_imax; i++){
            for(int k=n_kmin; k<n_kmax; k++){
                sum += conj(wf_buf[n][i][k])*(-1.0*_i[i]) * wf_buf[n][i][k]*_di*_dk;
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

    cdouble*** wf_buf = _wf->get_buf();
    for(int n=0;n<_nt_diag;n++){
        cdouble sum = 0.0; 
	#pragma omp parallel for schedule(dynamic) collapse(1)
        for(int i=n_imin; i<n_imax; i++){
            for(int k=n_kmin; k<n_kmax; k++){
                sum += conj(wf_buf[n][i][k])*(-1.0*_k[k]) * wf_buf[n][i][k]*_di*_dk;
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

    cdouble*** wf_buf = _wf->get_buf();
    for(int n=0;n<_nt_diag;n++){
        cdouble sum = 0.0; 
	#pragma omp parallel for schedule(dynamic) collapse(1)
        for(int i=n_imin; i<n_imax; i++){
            for(int k=n_kmin; k<n_kmax; k++){
                sum += conj(wf_buf[n][i][k])* wf_buf[n][i][k]*_di*_dk;
            }
        }
        _data[idx + n] = sum;
    }
}

void ProbeXZ::_pop_0(const int idx){
}

void ProbeXZ::_dens(const int idx){}
void ProbeXZ::_wf_snap(const int idx){}

ProbeXZ::~ProbeXZ(){
    delete _data;
}
