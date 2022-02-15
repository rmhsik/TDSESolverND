#include <iostream>
#include <tuple>
#include <string>
#include <cstdint>
#include <cstring>
#include <omp.h>
#include "tdsesolver.h"
#include "debug.h"
#include "utils.h"

void TDSESolver::_geom_XYZ(){
    std::tie(_i, _di) = linspace<double>(_param->imin, _param->imax, _param->ni);
    std::tie(_j, _dj) = linspace<double>(_param->jmin, _param->jmax, _param->nj);
    std::tie(_k, _dk) = linspace<double>(_param->kmin, _param->kmax, _param->nk);

    _propagate = &TDSESolver::_propagate_XYZ;
    _ipropagate = &TDSESolver::_ipropagate_XYZ;
}

void TDSESolver::_fields_XYZ(){
    std::string path;

    Afield_i = new Field(_param->E0i, _param->w0Ei, _param->phiEi, _param->env, _param->tmax_ev, _t, _param->nt);
    path = "results/Efield_i.dat";
    write_array(Afield_i->get(),_param->nt,path);
    Afield_i->calc_pot();
    path = "results/Afield_i.dat";
    write_array(Afield_i->get(), _param->nt, path);
    Bfield_i = new Field(_param->B0i, _param->w0Bi, _param->phiEi, _param->env, _param->tmax_ev, _t, _param->nt);
    path = "results/Bfield_i.dat";
    write_array(Bfield_i->get(), _param->nt, path);

    Afield_j = new Field(_param->E0j, _param->w0Ej, _param->phiEj, _param->env, _param->tmax_ev, _t, _param->nt);
    path = "results/Efield_j.dat";
    write_array(Afield_j->get(), _param->nt, path);
    Afield_j->calc_pot();
    path = "results/Afield_j.dat";
    write_array(Afield_j->get(), _param->nt, path);
    Bfield_j = new Field(_param->B0j, _param->w0Bj, _param->phiBj, _param->env, _param->tmax_ev, _t, _param->nt);
    path = "results/Bfield_j.dat";
    write_array(Bfield_j->get(), _param->nt, path);

    Afield_k = new Field(_param->E0k, _param->w0Ek, _param->phiEk, _param->env, _param->tmax_ev, _t, _param->nt);
    path = "results/Efield_k.dat";
    write_array(Afield_k->get(),_param->nt, path);
    Afield_k->calc_pot();
    path = "results/Afield_k.dat";
    write_array(Afield_k->get(),_param->nt, path);
    Bfield_k =new Field(_param->B0k, _param->w0Bk, _param->phiBk, _param->env, _param->tmax_ev, _t, _param->nt);
    path = "results/Bfield_k.dat";
    write_array(Bfield_k->get(),_param->nt, path);
}

void TDSESolver::_masks_XYZ(){
    double ib = _param->imax/10.0;
    double jb = _param->jmax/10.0;
    double kb = _param->kmax/10.0;
    double gamma = 1.0;

    for(int i=0; i<_param->ni; i++){
        if(_i[i] < _i[0]+ib){
            _imask[i] = pow(cos(M_PI*(_i[i]-(_i[0]+ib))*gamma/(2.0*ib)),1.0/8.0);
            if(_imask[i].real() < 0.0)
                _imask[i] = 0.0;
        }
        else if(_i[i]>(_i[_param->ni-1] - ib)){
            _imask[i] = pow(cos(M_PI*(_i[i] - (_i[_param->ni - 1]-ib))*gamma/(2.0*ib)),1.0/8.0);
        }
        else{
            _imask[i] = 1.0;
        }
    }

    for(int j=0; j<_param->nj; j++){
        if(_j[j] < _j[0]+jb){
            _jmask[j] = pow(cos(M_PI*(_j[j]-(_j[0]+jb))*gamma/(2.0*jb)),1.0/8.0);
            if(_jmask[j].real() < 0.0)
                _jmask[j] = 0.0;
        }
        else if(_j[j]>(_j[_param->nj-1] - jb)){
            _jmask[j] = pow(cos(M_PI*(_j[j] - (_j[_param->nj - 1]-jb))*gamma/(2.0*jb)),1.0/8.0);
        }
        else{
            _jmask[j] = 1.0;
        }
    }

    for(int k=0; k<_param->nk; k++){
        if(_k[k] < _k[0]+kb){
            _kmask[k] = pow(cos(M_PI*(_k[k]-(_k[0]+kb))*gamma/(2.0*kb)),1.0/8.0);
            if(_kmask[k].real() < 0.0)
                _kmask[k] = 0.0;
        }
        else if(_k[k]>(_k[_param->nk-1] - kb)){
            _kmask[k] = pow(cos(M_PI*(_k[k] - (_k[_param->nk - 1]-kb))*gamma/(2.0*kb)),1.0/8.0);
        }
        else{
            _kmask[k] = 1.0;
        }
    }
}

void TDSESolver::_ipropagate_XYZ(){
    cdouble ener = 0.0;
    cdouble norm;
    cdouble **psi_i_row, **psi_j_row, **psi_k_row;
    const int ni = _param->ni;
    const int nj = _param->nj;
    const int nk = _param->nk;
    psi_i_row = alloc2d<cdouble>(_param->n_threads,ni);
    psi_j_row = alloc2d<cdouble>(_param->n_threads,nj);
    psi_k_row = alloc2d<cdouble>(_param->n_threads,nk);

    for(int n=0; n<_param->nt_ITP;n++){
        double tstart, tend;
        tstart = omp_get_wtime();
        #pragma omp parallel for collapse(1) schedule(dynamic)
        for(int j=0;j<nj;j++){
            for(int k=0;k<nk;k++){
                int id = omp_get_thread_num();
                _wf->get_i_row(psi_i_row[id],j,k);
                (_ham->*(_ham->step_i))(psi_i_row[id],j,k,0,1,id);
                _wf->set_i_row(psi_i_row[id],j,k);
            }
        }

        #pragma omp parallel for collapse(1) schedule(dynamic)
        for(int i=0;i<ni;i++){
            for(int k=0;k<nk;k++){
                int id = omp_get_thread_num();
                _wf->get_j_row(psi_j_row[id],i,k);
                (_ham->*(_ham->step_j))(psi_j_row[id],i,k,0,1,id);
                _wf->set_j_row(psi_j_row[id],i,k);
            }
        }

        #pragma omp parallel for collapse(1) schedule(dynamic)
        for(int i=0;i<ni;i++){
            for(int j=0;j<nj;j++){
                int id = omp_get_thread_num();
                _wf->get_k_row(psi_k_row[id],i,j);
                (_ham->*(_ham->step_k))(psi_k_row[id],i,j,0,1,id);
                _wf->set_k_row(psi_k_row[id],i,j);
            }
        }
        norm = _wf->norm();
        (*_wf) /= norm;
        tend = omp_get_wtime();
        std::cout<<"n: "<<n<<" timestep: "<<tend-tstart<<"\n";
        if(n%5==0){
            ener = (_ham->*(_ham->ener))(_wf->get());
            std::cout<<"Norm: "<< norm<<" Ener: "<<ener<<"\n";
        }

    }
    std::cout<<"Ener: "<<ener<<"\n";
    free2d(&psi_i_row,_param->n_threads,ni);
    free2d(&psi_j_row,_param->n_threads,nj);
    free2d(&psi_k_row,_param->n_threads,nk);
}

void TDSESolver::_propagate_XYZ(){
    std::cout<<"Someday this might be working :("<<std::endl;
}
