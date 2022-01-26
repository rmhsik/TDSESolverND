#include <iostream>
#include <tuple>
#include <string>
#include <cstdint>
#include <cstring>
#include <omp.h>
#include "tdsesolver.h"
#include "debug.h"
#include "utils.h"

void TDSESolver::_geom_RZ(){
    std::tie(_i,_di) = linspace<double>(0.0,_param->imax,_param->ni);
    std::tie(_k,_dk) = linspace<double>(_param->kmin,_param->kmax, _param->nk);
    
    for(int i=0; i<_param->ni;i++){
        _i[i] += 0.5*_di;
    }

    _propagate = &TDSESolver::_propagate_RZ;
    _ipropagate = &TDSESolver::_ipropagate_RZ;
}

void TDSESolver::_fields_RZ(){
    std::string path;
    Afield_k = new Field(_param->E0k, _param->w0Ek, _param->phiEk, _param->env, _param->tmax_ev, _t, _param->nt);
    path = "results/Efield_k.dat";
    write_array(Afield_k->get(),_param->nt,path);
    Bfield_k = new Field(_param->B0k, _param->w0Bk, _param->phiBk, _param->env, _param->tmax_ev, _t, _param->nt);
    Afield_k->calc_pot();
    path = "results/Afield_k.dat";
    write_array(Afield_k->get(),_param->nt,path);
    path = "results/Bfield_k.dat";
    write_array(Bfield_k->get(),_param->nt,path);

    Afield_i = NULL;
    Bfield_i = NULL;
}

void TDSESolver::_masks_RZ(){
    double ib = _param->imax/10.0;
    double kb = _param->kmax/10.0;
    double gamma = 1.0;
    std::string path;

    for(int i=0; i<_param->ni;i++){
        if (_i[i]>(_i[_param->ni-1]-ib)){
            _imask[i] = pow(cos(M_PI*(_i[i] - (_i[_param->ni-1]-ib))*gamma/(2*ib)),1.0/8.0);
        }
        else{
            _imask[i] = 1.0;
        }
    }

    for(int i=0; i<_param->nk;i++){
        if(_k[i]<_k[0]+kb){
            _kmask[i] = pow(cos(M_PI*(_k[i] - (_k[0]+kb))*gamma/(2*kb)),1.0/8.0);
        }
        else if (_k[i]>(_k[_param->nk-1]-kb)){
            _kmask[i] = pow(cos(M_PI*(_k[i] - (_k[_param->nk-1]-kb))*gamma/(2*kb)),1.0/8.0);
        }
        else{
            _kmask[i] = 1.0;
        }
    }
    path = "results/imask.dat";
    write_array(_imask,_param->ni,path);
    path = "results/kmask.dat";
    write_array(_kmask,_param->nk,path);
}


void TDSESolver::_ipropagate_RZ(){
    cdouble ener = 0.0;
    cdouble norm;
    cdouble *psi_col, *psi_row;
    const int ni = _param->ni;
    const int nk = _param->nk;
    psi_col = new cdouble [_param->nk*_param->n_threads];
    psi_row = new cdouble [_param->ni*_param->n_threads];

    for(int j=0; j<_param->nt_ITP;j++){
        #pragma omp parallel for schedule(dynamic)
        for(int i=0;i<ni;i++){
            int id_thread = omp_get_thread_num();
            for(int k=0;k<nk;k++)
                psi_col[id_thread*nk + k] = _wf->get()[i*nk+k];
            (_ham->*(_ham->step_k))(&psi_col[id_thread*nk],0.0,0.0,i,1,id_thread);
            _wf->set_col(&psi_col[id_thread*nk],i);
        }
        
        #pragma omp parallel for schedule(dynamic)
        for(int k=0;k<nk;k++){
            int id_thread = omp_get_thread_num();
            for(int i=0;i<ni;i++)
                psi_row[id_thread*ni + i] = _wf->get()[i*nk+k];
            (_ham->*(_ham->step_i))(&psi_row[id_thread*ni],0.0,0.0,k,1,id_thread);
            _wf->set_row(&psi_row[id_thread*ni],k);
        }
        
        #pragma omp parallel for schedule(dynamic)
        for(int i=0;i<ni;i++){
            int id_thread = omp_get_thread_num();
            for(int k=0;k<nk;k++)
                psi_col[id_thread*nk + k] = _wf->get()[i*nk+k];
            (_ham->*(_ham->step_k))(&psi_col[id_thread*nk],0.0,0.0,i,1,id_thread);
            _wf->set_col(&psi_col[id_thread*nk],i);
        }

        norm = _wf->norm();
        (*_wf) /= norm;
        if(j%50==0){
            ener = (_ham->*(_ham->ener))(_wf->get());
            std::cout<<"Norm: "<< norm<<" Ener: "<<ener<<"\n";
        }
    }
    std::cout<<"Ener: "<<ener<<"\n";
    delete psi_col;
    delete psi_row;
}

void TDSESolver::_propagate_RZ(){
    cdouble norm, ener;
    cdouble *acc_i_vec, *acc_k_vec, *dip_vec, *pop_vec;
    cdouble *psi_col, *psi_row;
    int idx;
    const int ni = _param->ni;
    const int nk = _param->nk;
    std::string path;

    acc_i_vec = new cdouble [_param->nt];
    acc_k_vec = new cdouble [_param->nt];
    dip_vec = new cdouble [_param->nt];
    pop_vec = new cdouble [_param->nt];

    psi_col = new cdouble [_param->nk*_param->n_threads];
    psi_row = new cdouble [_param->ni*_param->n_threads];

    //wf_ptr = _wf.get_buf();
    _wf->set_to_buf(0); 
    for(int j=0; j<_param->nt;j++){
        #pragma omp parallel for schedule(dynamic)
        for(int i=0;i<ni;i++){
            cdouble *wf_ptr;
            wf_ptr = _wf->get_buf();
            int id_thread = omp_get_thread_num();
            for(int k=0;k<nk;k++)
                psi_col[id_thread*nk + k] = wf_ptr[j%_param->nt_diag*ni*nk+i*nk+k];
            (_ham->*(_ham->step_k))(&psi_col[id_thread*nk],(*Afield_k)[j],(*Bfield_k)[j],i,0,id_thread);
            _wf->set_col_buf_mask(&psi_col[id_thread*nk],_kmask, i,(j+1)%_param->nt_diag);
        }
        
        #pragma omp parallel for schedule(dynamic)
        for(int k=0;k<nk;k++){
            cdouble *wf_ptr;
            wf_ptr = _wf->get_buf();
            int id_thread = omp_get_thread_num();
            for(int i=0;i<ni;i++)
                psi_row[id_thread*ni + i] = wf_ptr[(j+1)%_param->nt_diag*ni*nk+i*nk+k];
            (_ham->*(_ham->step_i))(&psi_row[id_thread*ni],(*Afield_k)[j],(*Bfield_k)[j],k,0,id_thread);
            _wf->set_row_buf_mask(&psi_row[id_thread*ni],_imask,k,(j+1)%_param->nt_diag);
        }
        
        #pragma omp parallel for schedule(dynamic)
        for(int i=0;i<ni;i++){
            cdouble *wf_ptr;
            wf_ptr = _wf->get_buf();
            int id_thread = omp_get_thread_num();
            for(int k=0;k<nk;k++)
                psi_col[id_thread*nk + k] = wf_ptr[(j+1)%_param->nt_diag*ni*nk+i*nk+k];
            (_ham->*(_ham->step_k))(&psi_col[id_thread*nk],(*Afield_k)[j],(*Bfield_k)[j],i,0,id_thread);
            _wf->set_col_buf_mask(&psi_col[id_thread*nk],_kmask,i,(j+1)%_param->nt_diag);
        }
        //_wf.apply_mask_buf_RZ(_imask,_kmask,(j+1)%_param->nt_diag);
        //_wf.set_to_buf(j%_param->nt_diag);
        //acc_vec[j] = _wf.acc();
        //dip_vec[j] = _wf.dipole();
        if ((j+2)%_param->nt_diag==0 && j<(_param->nt-_param->nt%_param->nt_diag)){ 
            idx = j - _param->nt_diag + 2;
            _wf->dip_k_buf();
            std::memcpy(&dip_vec[idx], _wf->get_diag_buf(), _param->nt_diag*sizeof(cdouble));
            _wf->acc_i_buf();
            std::memcpy(&acc_i_vec[idx], _wf->get_diag_buf(), _param->nt_diag*sizeof(cdouble));
            _wf->acc_k_buf();
            std::memcpy(&acc_k_vec[idx], _wf->get_diag_buf(), _param->nt_diag*sizeof(cdouble));
            _wf->pop_buf(_param->pop_imin, _param->pop_imax, _param->pop_kmin, _param->pop_kmax); 
            std::memcpy(&pop_vec[idx], _wf->get_diag_buf(), _param->nt_diag*sizeof(cdouble));
            //norm = _wf.norm();
            //ener = (_ham->*(_ham->ener))(_wf.get());
            //std::cout<<j<<" Norm: "<< norm<<"\n";
        }
    }
    // Get last batch of diagnostics from last idx up to _paran.nt
    _wf->dip_k_buf();
    std::memcpy(&dip_vec[idx+1],_wf->get_diag_buf(),(_param->nt%_param->nt_diag-1)*sizeof(cdouble));
    _wf->acc_i_buf();
    std::memcpy(&acc_i_vec[idx+1],_wf->get_diag_buf(),(_param->nt%_param->nt_diag-1)*sizeof(cdouble));
    _wf->acc_k_buf();
    std::memcpy(&acc_k_vec[idx+1],_wf->get_diag_buf(),(_param->nt%_param->nt_diag-1)*sizeof(cdouble));
    _wf->pop_buf(_param->pop_imin, _param->pop_imax, _param->pop_kmin, _param->pop_kmax);
    std::memcpy(&pop_vec[idx], _wf->get_diag_buf(), (_param->nt%_param->nt_diag-1)*sizeof(cdouble));
    cdouble valaccmask;
    for(int j=0; j<_param->nt; j++){
        valaccmask = _accmask[j];
        acc_i_vec[j] *= valaccmask;
        acc_k_vec[j] *= valaccmask;
        dip_vec[j] *= valaccmask;
    }

    write_array(acc_i_vec, _param->nt, _param->acc_i_path);
    write_array(acc_k_vec, _param->nt, _param->acc_k_path);
    write_array(dip_vec, _param->nt, _param->dip_path);
    write_array(pop_vec, _param->nt, _param->pop_path);
    delete acc_i_vec;
    delete acc_k_vec;
    delete dip_vec;
    delete pop_vec;
}


