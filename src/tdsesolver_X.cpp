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
    std::tie(_j,_dj) = linspace<double>(_param->jmin,_param->jmax,_param->nj);
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
    _jmask[0] = 1.0;
    path = "results/imask.dat";
    write_array(_imask,_param->ni,path);
    path = "results/jmask.dat";
    write_array(_jmask,_param->nj,path);
    path = "results/kmask.dat";
    write_array(_kmask,_param->nk,path);

}

void TDSESolver::_ipropagate_X(){
    debug3("[TDSESolver->ipropagate] Start imaginary propagation...");
    cdouble norm;    
    std::cout<<"Norm: "<<_wf->norm()<<"\n";   
    for(int i=0; i<_param->nt_ITP; i++){
        cdouble *psi_row;
        psi_row = _wf->i_row(0,0);

        (_ham->*(_ham->step_i))(psi_row,0,0,0,1,0);
        _wf->set_i_row(psi_row,0,0);
        norm = _wf->norm();
        (*_wf)/=norm;
        if(i%100 ==0){
            std::cout<<"Norm:"<<norm<<" "<<(_ham->*(_ham->ener))(_wf->get())<<"\n";
        }
    }

    std::cout<<"Norm: "<<_wf->norm()<<"\n";   
    std::string path = "results/itp_psi2.dat";
    
    debug3("[TDSESolver->ipropagate] End imaginary propagation");


}

void TDSESolver::_propagate_X(){
    debug3("[TDSESolver->propagate] Start propagate...");
    cdouble *psi_row;
    cdouble *norm_vec;
    int norm_vec_size = (int)(_param->nt/100);
    int norm_vec_idx;
    int idx;
    _wf->set_to_buf(0);
    for(int i=0; i<_param->nt; i++){
        psi_row = _wf->i_row(0,0);
        (_ham->*(_ham->step_i))(psi_row,0,0,i,0,0);
        _wf->set_i_row(psi_row,0,0);
        _wf->apply_mask(_imask,_jmask,_kmask);
        _wf->set_to_buf(i%_param->nt_diag);

        if((i+1)%_param->nt_diag==0 && i<(_param->nt - _param->nt%_param->nt_diag)){
            idx = i - _param->nt_diag + 1;
            _diag->run_diagnostics(idx);
        }
    }
    // Get las batch if diagnostics from lat idx up to _param.nt
    _diag->write_diagnostics();    
    debug3("[TDSESolver->propagate] End propagate");
}

