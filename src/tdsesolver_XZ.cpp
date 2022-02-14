#include <iostream>
#include <tuple>
#include <string>
#include <cstdint>
#include <cstring>
#include <omp.h>
#include "tdsesolver.h"
#include "debug.h"
#include "utils.h"
#include "diagnostics.h"

void TDSESolver::_geom_XZ(){
    std::tie(_i,_di) = linspace<double>(_param->imin,_param->imax,_param->ni);
    std::tie(_j,_dj) = linspace<double>(_param->jmin,_param->jmax,_param->nj);
    std::tie(_k,_dk) = linspace<double>(_param->kmin,_param->kmax,_param->nk);
    _propagate = &TDSESolver::_propagate_XZ;
    _ipropagate = &TDSESolver::_ipropagate_XZ;
}

void TDSESolver::_fields_XZ(){
    std::string path;
    Afield_i = new Field(_param->E0i, _param->w0Ei, _param->phiEi, _param->env, _param->tmax_ev, _t, _param->nt);
    Afield_k = new Field(_param->E0k, _param->w0Ek, _param->phiEk, _param->env, _param->tmax_ev, _t, _param->nt);
    Bfield_i = new Field(_param->B0i, _param->w0Bi, _param->phiBi, _param->env, _param->tmax_ev, _t, _param->nt);
    Bfield_k = new Field(_param->B0k, _param->w0Bk, _param->phiBk, _param->env, _param->tmax_ev, _t, _param->nt);
    Afield_i->calc_pot();
    Afield_k->calc_pot();
    path = "results/Afield_i.dat";
    write_array(Afield_i->get(),_param->nt,path);
    path = "results/Afield_k.dat";
    write_array(Afield_k->get(),_param->nt,path);
    path = "results/Bfield_i.dat";
    write_array(Bfield_i->get(),_param->nt,path);
    path = "results/Bfield_k.dat";
    write_array(Bfield_k->get(),_param->nt,path);
}

void TDSESolver::_masks_XZ(){
    double ib = _param->imax/10.0;
    double kb = _param->kmax/10.0;
    double gamma = 1.0;
    std::string path;

    for(int i=0; i<_param->ni;i++){
        if(_i[i]<_i[0]+ib){
            _imask[i] = pow(cos(M_PI*(_i[i] - (_i[0]+ib))*gamma/(2*ib)),1.0/8.0);
            if(_imask[i].real()<0.0)
                _imask[i] =0.0;
        }
        else if (_i[i]>(_i[_param->ni-1]-ib)){
            _imask[i] = pow(cos(M_PI*(_i[i] - (_i[_param->ni-1]-ib))*gamma/(2*ib)),1.0/8.0);
            if(_imask[i].real()<0.0)
                _imask[i] =0.0;

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
    _jmask[0] = 1.0;
    path = "results/imask.dat";
    write_array(_imask,_param->ni,path);
    path = "results/jmask.dat";
    write_array(_jmask,_param->nj,path);
    path = "results/kmask.dat";
    write_array(_kmask,_param->nk,path);
}

void TDSESolver::_ipropagate_XZ(){
    cdouble ener = 0.0;
    cdouble norm;
    cdouble **psi_i_row, **psi_k_row;
    const int ni = _param->ni;
    const int nk = _param->nk;
    psi_i_row = alloc2d<cdouble>(_param->n_threads, ni);
    psi_k_row = alloc2d<cdouble>(_param->n_threads, nk);

    for(int n=0; n<_param->nt_ITP;n++){
        #pragma omp parallel for schedule(dynamic)
        for(int i=0;i<ni;i++){
            int id_thread = omp_get_thread_num();
            _wf->get_k_row(psi_k_row[id_thread],i,0);
            (_ham->*(_ham->step_k))(psi_k_row[id_thread],i,0,0,1,id_thread);
            _wf->set_k_row(psi_k_row[id_thread],i,0);
        }
        
        #pragma omp parallel for schedule(dynamic)
        for(int k=0;k<nk;k++){
            int id_thread = omp_get_thread_num();
            _wf->get_i_row(psi_i_row[id_thread],0,k);
            (_ham->*(_ham->step_i))(psi_i_row[id_thread],0,k,0,1,id_thread);
            _wf->set_i_row(psi_i_row[id_thread],0,k);
        }
        
        norm = _wf->norm();
        (*_wf) /= norm;
        if(n%50==0){
            ener = (_ham->*(_ham->ener))(_wf->get());
            std::cout<<"Norm: "<< norm<<" Ener: "<<ener<<"\n";
        }
    }
    std::cout<<"Ener: "<<ener<<"\n";
    free2d(&psi_i_row,_param->n_threads,ni);
    free2d(&psi_k_row,_param->n_threads,nk);
}


void TDSESolver::_propagate_XZ(){
    cdouble norm, ener;
    cdouble **psi_i_row, **psi_k_row;
    int idx;
    const int ni = _param->ni;
    const int nk = _param->nk;

    psi_i_row = alloc2d<cdouble>(_param->n_threads, ni);
    psi_k_row = alloc2d<cdouble>( _param->n_threads, nk);

    //wf_ptr = _wf.get_buf();
    _wf->set_to_buf(0); 
    for(int n=0; n<_param->nt;n++){
        #pragma omp parallel for schedule(dynamic)
        for(int i=0;i<ni;i++){
            int id_thread = omp_get_thread_num();
            _wf->get_k_row_buf(psi_k_row[id_thread],i,0,n%_param->nt_diag);
            (_ham->*(_ham->step_k))(psi_k_row[id_thread],i,0,n,0,id_thread);
            _wf->set_k_row_buf_mask(psi_k_row[id_thread],_kmask,i,0,(n+1)%_param->nt_diag);
        }
        
        #pragma omp parallel for schedule(dynamic)
        for(int k=0;k<nk;k++){
            int id_thread = omp_get_thread_num();
            _wf->get_i_row_buf(psi_i_row[id_thread],0,k,(n+1)%_param->nt_diag);
            (_ham->*(_ham->step_i))(psi_i_row[id_thread],0,k,n,0,id_thread);
            _wf->set_i_row_buf_mask(psi_i_row[id_thread],_imask,0,k,(n+1)%_param->nt_diag);
        }

        if ((n+2)%_param->nt_diag==0 && n<(_param->nt-_param->nt%_param->nt_diag)){ 
             idx = n - _param->nt_diag + 2;
             _diag->run_diagnostics(idx);
        }
    }
    //TODO: Get last batch of diagnostics from last idx up to _paran.nt
    _diag->write_diagnostics();
    free2d(&psi_i_row,_param->n_threads,ni);
    free2d(&psi_k_row,_param->n_threads,nk);
    
}
