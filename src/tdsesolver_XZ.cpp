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
    std::tie(_k,_dk) = linspace<double>(_param->kmin,_param->kmax,_param->nk);
    _propagate = &TDSESolver::_propagate_XZ;
    _ipropagate = &TDSESolver::_ipropagate_XZ;
}

void TDSESolver::_fields_XZ(){
    std::string path;
    //Afield_i = new Field(_param->E0i, _param->w0Ei, _param->phiEi, _param->env, _param->tmax_ev, _t, _param->nt);
    Afield_k = new Field(_param->E0k, _param->w0Ek, _param->kk, _param->phiEk, _param->env, _param->tmax_ev, _t, _param->nt);
    Afield_k->set_params(_param);
    Afield_k->set_geometry(_i,_k,_di,_dk);
    //Afield_i->calc_pot();
    //Afield_k->calc_pot();
    //path = "results/Afield_i.dat";
    //write_array(Afield_i->get(),_param->nt,path);
    //path = "results/Afield_k.dat";
    //write_array(Afield_k->get(),_param->nt,path);
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
    path = "results/imask.dat";
    write_array(_imask,_param->ni,path);
    path = "results/kmask.dat";
    write_array(_kmask,_param->nk,path);
}

void TDSESolver::_ipropagate_XZ(){

    if (_mpi_grid->rank == 0){
        cdouble ener = 0.0;
        cdouble norm;
        cdouble **psi_i_row, **psi_k_row;
        const int ni = _param->ni;
        const int nk = _param->nk;
        psi_i_row = alloc2d<cdouble>(_param->n_threads, ni);
        psi_k_row = alloc2d<cdouble>(_param->n_threads, nk);
        for(int n=0; n<_param->nt_ITP;n++){
            //#pragma omp parallel for schedule(dynamic)
            for(int i=0;i<ni;i++){
                int id_thread = 0;//omp_get_thread_num();
                _wf->get_k_row_mpi(psi_k_row[id_thread],i);
                (_ham->*(_ham->step_k))(psi_k_row[id_thread],i,0,1,id_thread);
                _wf->set_k_row(psi_k_row[id_thread],i);
            }
            
            //#pragma omp parallel for schedule(dynamic)
            for(int k=0;k<nk;k++){
                int id_thread = 0;//omp_get_thread_num();
                _wf->get_i_row(psi_i_row[id_thread],k);
                (_ham->*(_ham->step_i))(psi_i_row[id_thread],k,0,1,id_thread);
                _wf->set_i_row(psi_i_row[id_thread],k);
            }
            
            norm = _wf->norm();
            (*_wf) /= norm;
            if(n%50==0){
                ener = (_ham->*(_ham->ener))(_wf->get());
                std::cout<<"Norm: "<< norm<<" Ener: "<<ener<<std::endl;
            }
        }
        std::cout<<"Ener: "<<ener<<std::endl;
        free2d(&psi_i_row,_param->n_threads,ni);
        free2d(&psi_k_row,_param->n_threads,nk);
    }
}


void TDSESolver::_propagate_XZ(){
    cdouble norm, ener;
    cdouble **psi_i_row, **psi_k_row;
    int idx;
    const int ni = _param->ni;
    const int nk = _param->nk;
    double *vecpot;
    vecpot = new double[_param->nt];

    psi_i_row = alloc2d<cdouble>(_param->n_threads, ni);
    psi_k_row = alloc2d<cdouble>( _param->n_threads, nk);

    //wf_ptr = _wf.get_buf();
    _wf->set_to_buf(0); 
    for(int n=0; n<_param->nt;n++){
        Afield_k->calc_field(_t[n]);
        Afield_k->calc_pot();
        //#pragma omp parallel for schedule(dynamic)
        for(int i=0;i<ni;i++){
            int id_thread = omp_get_thread_num();
            _wf->get_k_row_buf(psi_k_row[id_thread],i,n%_param->nt_diag);
            (_ham->*(_ham->step_k))(psi_k_row[id_thread],i,n,0,id_thread);
            _wf->set_k_row_buf_mask(psi_k_row[id_thread],_kmask,i,(n+1)%_param->nt_diag);
        }
        
        //#pragma omp parallel for schedule(dynamic)
        for(int k=0;k<nk;k++){
            int id_thread = omp_get_thread_num();
            _wf->get_i_row_buf(psi_i_row[id_thread],k,(n+1)%_param->nt_diag);
            (_ham->*(_ham->step_i))(psi_i_row[id_thread],k,n,0,id_thread);
            _wf->set_i_row_buf_mask(psi_i_row[id_thread],_imask,k,(n+1)%_param->nt_diag);
        }
        
        if ((n+2)%_param->nt_diag==0 && n<(_param->nt-_param->nt%_param->nt_diag)){ 
             idx = n - _param->nt_diag + 2;
             _diag->run_diagnostics(idx);
             norm = _wf->norm_buf(_param->nt_diag-1);
             std::cout<<"ti: "<< n<<" adield: "<< norm<<std::endl;
        }
        vecpot[n] = Afield_k->get(_param->ni/2,_param->nk/2); 
    }
    //TODO: Get last batch of diagnostics from last idx up to _paran.nt
    //
    _diag->write_diagnostics();
    free2d(&psi_i_row,_param->n_threads,ni);
    free2d(&psi_k_row,_param->n_threads,nk);
    write_array(vecpot,_param->nt,"results/vecpot.dat");
}
