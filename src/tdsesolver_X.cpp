#include <iostream>
#include <tuple>
#include <string>
#include <cstdint>
#include <cstring>
#include <omp.h>
#include "tdsesolver.h"
#include "debug.h"
#include "utils.h"

void TDSESolver::_geom_X(){
    std::tie(_i,_di) = linspace<double>(_param->imin,_param->imax,_param->ni);
    std::tie(_k,_dk) = linspace<double>(_param->kmin,_param->kmax,_param->nk);
    _propagate = &TDSESolver::_propagate_X;
    _ipropagate = &TDSESolver::_ipropagate_X;

}

void TDSESolver::_fields_X(){
    std::string path;
    Afield_i = new Field(_param->E0i, _param->w0Ei, _param->phiEi, _param->env, _param->tmax_ev,_t, _param->nt);
    Afield_i->calc_pot();
    path = "results/Afield_i.dat";
    write_array(Afield_i->get(),_param->nt,path);
    Afield_k = NULL;
    Bfield_i = NULL;
    Bfield_k = NULL;
}

void TDSESolver::_masks_X(){
    double ib = _param->imax/10.0;
    double gamma = 1.0;
    std::string path;

    for(int i=0; i<_param->ni;i++){
        if(_i[i]<_i[0]+ib){
            _imask[i] = pow(cos(M_PI*(_i[i] - (_i[0]+ib))*gamma/(2*ib)),1.0/8.0);
        }
        else if (_i[i]>(_i[_param->ni-1]-ib)){
            _imask[i] = pow(cos(M_PI*(_i[i] - (_i[_param->ni-1]-ib))*gamma/(2*ib)),1.0/8.0);
        }
        else{
            _imask[i] = 1.0;
        }
    }
    _kmask[0] = 1.0;
    path = "results/imask.dat";
    write_array(_imask,_param->ni,path);
    path = "results/kmask.dat";
    write_array(_kmask,_param->nk,path);

}

void TDSESolver::_ipropagate_X(){
    debug3("[TDSESolver->ipropagate] Start imaginary propagation...");
    cdouble norm;    
    for(int i=0; i<_param->nt_ITP; i++){
        cdouble *psi_row;
        psi_row = _wf->row(0);

        (_ham->*(_ham->step_i))(psi_row,0.0,0.0,0,1,0);
        //(_ham->.*(_ham->step_i)(psi_row,0.0,0.0,1);
        _wf->set_row(psi_row,0);
        norm = _wf->norm();
        (*_wf)/=norm;
        if(i%100 ==0){
            std::cout<<"Norm:"<<norm<<" "<<(_ham->*(_ham->ener))(_wf->get())<<"\n";
        }
    }

    std::cout<<"Norm: "<<_wf->norm()<<"\n";   
    std::string path = "results/itp_psi2.dat";
    _wf->save_wf2(path);
    
    debug3("[TDSESolver->ipropagate] End imaginary propagation");


}

void TDSESolver::_propagate_X(){
    debug3("[TDSESolver->propagate] Start propagate...");
    cdouble *psi_row, *dip_vec,*acc_vec, *pop_vec;
    cdouble *norm_vec;
    int norm_vec_size = (int)(_param->nt/100);
    int norm_vec_idx;
    int idx;
    acc_vec = new cdouble[_param->nt];
    dip_vec = new cdouble[_param->nt];
    pop_vec = new cdouble[_param->nt];
    //norm_vec = new cdouble[norm_vec_size];
    for(int i=0; i<_param->nt; i++){
        psi_row = _wf->row(0);
        (_ham->*(_ham->step_i))(psi_row,(*Afield_i)[i],0.0,0,0,0);
        _wf->set_row(psi_row,0);
        _wf->apply_mask(_imask,_kmask);
        
        _wf->set_to_buf(i%_param->nt_diag);
        //acc_vec[i] = _wf.acc();
        if ((i+2)%_param->nt_diag==0 && i<(_param->nt-_param->nt%_param->nt_diag)){ 
            idx = i - _param->nt_diag + 2;
            _wf->acc_k_buf();
            std::memcpy(&acc_vec[idx],_wf->get_diag_buf(),_param->nt_diag*sizeof(cdouble));
            _wf->dip_i_buf();
            std::memcpy(&dip_vec[idx],_wf->get_diag_buf(),_param->nt_diag*sizeof(cdouble));
            _wf->pop_buf(_param->pop_imin, _param->pop_imax, _param->pop_kmin, _param->pop_kmax);
            std::memcpy(&pop_vec[idx],_wf->get_diag_buf(),_param->nt_diag*sizeof(cdouble));
            //norm_vec_idx = (int)(i/100);
            //norm_vec[norm_vec_idx] = _wf.norm();
            //std::cout<<"Ener: "<<(_ham->*(_ham->ener))(psi_row)<<"\n";
        }
    }
    // Get las batch if diagnostics from lat idx up to _param.nt
    _wf->acc_k_buf();
    std::memcpy(&acc_vec[idx+1],_wf->get_diag_buf(),(_param->nt%_param->nt_diag-1)*sizeof(cdouble));
    _wf->dip_i_buf();
    std::memcpy(&dip_vec[idx+1],_wf->get_diag_buf(),(_param->nt%_param->nt_diag-1)*sizeof(cdouble));
    _wf->pop_buf(_param->pop_imin, _param->pop_imax, _param->pop_kmin, _param->pop_kmax);
    std::memcpy(&pop_vec[idx+1],_wf->get_diag_buf(),(_param->nt%_param->nt_diag-1)*sizeof(cdouble));

    
    for(int i=0; i<_param->nt;i++){
        acc_vec[i] *= _accmask[i];
        dip_vec[i] *= _accmask[i];
    }
    write_array(acc_vec,_param->nt,_param->acc_k_path);
    write_array(dip_vec,_param->nt,_param->dip_i_path);
    write_array(pop_vec,_param->nt,_param->pop_path);
   
    delete acc_vec;
    delete dip_vec;
    delete pop_vec;
    debug3("[TDSESolver->propagate] End propagate");
}

