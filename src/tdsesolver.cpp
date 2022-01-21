#include <iostream>
#include <tuple>
#include <string>
#include <cstdint>
#include <cstring>
#include <omp.h>
#include "tdsesolver.h"
#include "debug.h"
#include "utils.h"

TDSESolver::TDSESolver(){
}

TDSESolver::TDSESolver(Parameters *param){
    std::string path;
    _param = param;
    setup_time();

    omp_set_num_threads(_param->n_threads);
    setup_geometry();
    setup_fields();
    setup_ham();
    setup_wf();
    setup_masks();
    }

void TDSESolver::setup_time(){
    std::string path;
    std::tie(_t,_dt) = linspace<double>(0.0,_param->tmax_sim, _param->nt);
    path = "results/time.dat";
    write_array(_t,_param->nt,path);
}

void TDSESolver::setup_geometry(){
    std::string path;
    switch(_param->geometry){
        case X:
            std::tie(_i,_di) = linspace<double>(_param->imin,_param->imax,_param->ni);
            std::tie(_k,_dk) = linspace<double>(_param->kmin,_param->kmax,_param->nk);
	        _propagate = &TDSESolver::propagate_X;
	        _ipropagate = &TDSESolver::ipropagate_X;
            break;
        case XZ:
            std::tie(_i,_di) = linspace<double>(_param->imin,_param->imax,_param->ni);
            std::tie(_k,_dk) = linspace<double>(_param->kmin,_param->kmax,_param->nk);
            _propagate = &TDSESolver::propagate_XZ;
	        _ipropagate = &TDSESolver::ipropagate_XZ;

            break;
        case RZ:
            std::tie(_i,_di) = linspace<double>(0.0,_param->imax,_param->ni);
            std::tie(_k,_dk) = linspace<double>(_param->kmin,_param->kmax, _param->nk);
            
            for(int i=0; i<_param->ni;i++){
                _i[i] += 0.5*_di;
            }

	        _propagate = &TDSESolver::propagate_RZ;
	        _ipropagate = &TDSESolver::ipropagate_RZ;
            break;
    }
    path = "results/i.dat";
    write_array(_i,_param->ni,path);
    path = "results/k.dat";
    write_array(_k,_param->nk,path);
}

void TDSESolver::setup_fields(){
    std::string path;
    switch(_param->geometry){
        case X:
            Afield_i = Field(_param->E0i, _param->w0Ei, _param->phiEi, _param->env, _param->tmax_ev,_t, _param->nt);
            Afield_i.calc_pot();
            path = "results/Afield_i.dat";
            write_array(Afield_i.get(),_param->nt,path);
            break;

        case XZ:
            Afield_i = Field(_param->E0i, _param->w0Ei, _param->phiEi, _param->env, _param->tmax_ev, _t, _param->nt);
            Afield_k = Field(_param->E0k, _param->w0Ek, _param->phiEk, _param->env, _param->tmax_ev, _t, _param->nt);
            Bfield_i = Field(_param->B0i, _param->w0Bi, _param->phiBi, _param->env, _param->tmax_ev, _t, _param->nt);
            Bfield_k = Field(_param->B0k, _param->w0Bk, _param->phiBk, _param->env, _param->tmax_ev, _t, _param->nt);
            Afield_i.calc_pot();
            Afield_k.calc_pot();
            path = "results/Afield_i.dat";
            write_array(Afield_i.get(),_param->nt,path);
            path = "results/Afield_k.dat";
            write_array(Afield_k.get(),_param->nt,path);
            path = "results/Bfield_i.dat";
            write_array(Bfield_i.get(),_param->nt,path);
            path = "results/Bfield_k.dat";
            write_array(Bfield_k.get(),_param->nt,path);
            break;

        case RZ:
            Afield_k = Field(_param->E0k, _param->w0Ek, _param->phiEk, _param->env, _param->tmax_ev, _t, _param->nt);
            path = "results/Efield_k.dat";
            write_array(Afield_k.get(),_param->nt,path);
            Bfield_k = Field(_param->B0k, _param->w0Bk, _param->phiBk, _param->env, _param->tmax_ev, _t, _param->nt);
            Afield_k.calc_pot();
            path = "results/Afield_k.dat";
            write_array(Afield_k.get(),_param->nt,path);
            path = "results/Bfield_k.dat";
            write_array(Bfield_k.get(),_param->nt,path);
            break;
    }
}

void TDSESolver::setup_wf(){
    _wf = WF(_param);
    _wf.set_geometry(_i,_k,_di,_dk);
    _wf.set_dpotential(_ham.get_dpotential()); 
    switch(_param->init_wf){
    	case GAUS:
            _wf.gaussian(0.0,0.0,1.0);
	    break;
	case EXPO:
	    _wf.exponential(0.0,0.0,1.0);
	    break;
    } 
    cdouble norm = _wf.norm();
    //std::cout<<norm<<std::endl;
    _wf /= norm;
    norm = _wf.norm();
    //std::cout<<norm<<std::endl;

    std::string path = "results/init_psi2.dat";
    //write_array(_wf.get(),_param->ni*_param->nk,path);
}

void TDSESolver::setup_ham(){
    _ham = Hamiltonian(_param);
    _ham.set_geometry(_i,_k,_di, _dk);
    _ham.set_potential();
    _ham.set_dpotential();
}

void TDSESolver::setup_masks(){
    _accmask = new cdouble[_param->nt];
    _imask = new cdouble[_param->ni];
    _kmask = new cdouble[_param->nk];
    double tmax_mask = 0.0;
    double ib = _param->imax/10.0;
    double kb = _param->kmax/10.0;
    double gamma = 1.0;
    std::string path;
    
    if (_param->env==SIN2){
        tmax_mask = _param->tmax_ev;
    }
    else{
        tmax_mask = _param->tmax_ev + _param->period*2.0;
    }

    for(int i=0; i<_param->nt;i++){
        if(_t[i]<tmax_mask){
            _accmask[i] = 1.0;
        }
        else{
            _accmask[i] = exp(-pow(_t[i]-tmax_mask,2)/100.0);
        }
    }
    path = "results/accmask.dat";
    write_array(_accmask,_param->nt,path);

    switch(_param->geometry){
        case X:
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

            break;

        case XZ:
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

            break; 

        case RZ:
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

            break;
    }

}

void TDSESolver::ipropagate(){
    (this->*(this->_ipropagate))();
}

void TDSESolver::propagate(){
    (this->*(this->_propagate))();
}

void TDSESolver::ipropagate_X(){
    debug3("[TDSESolver->ipropagate] Start imaginary propagation...");
    cdouble norm;    
    for(int i=0; i<_param->nt_ITP; i++){
        cdouble *psi_row;
        psi_row = _wf.row(0);

        (_ham.*(_ham.step_i))(psi_row,0.0,0.0,0,1,0);
        //(_ham).*(_ham.step_i)(psi_row,0.0,0.0,1);
        _wf.set_row(psi_row,0);
        norm = _wf.norm();
        _wf/=norm;
        if(i%100 ==0){
            std::cout<<"Norm:"<<norm<<" "<<(_ham.*(_ham.ener))(_wf.get())<<"\n";
        }
    }

    std::cout<<"Norm: "<<_wf.norm()<<"\n";   
    std::string path = "results/itp_psi2.dat";
    _wf.save_wf2(path);
    
    debug3("[TDSESolver->ipropagate] End imaginary propagation");


}

void TDSESolver::propagate_X(){
    debug3("[TDSESolver->propagate] Start propagate...");
    cdouble *psi_row, *acc_vec;
    cdouble *norm_vec;
    int norm_vec_size = (int)(_param->nt/100);
    int norm_vec_idx;

    acc_vec = new cdouble[_param->nt];
    //norm_vec = new cdouble[norm_vec_size];
    for(int i=0; i<_param->nt; i++){
        psi_row = _wf.row(0);
        (_ham.*(_ham.step_i))(psi_row,Afield_i[i],0.0,0,0,0);
        _wf.set_row(psi_row,0);
        _wf.apply_mask(_imask,_kmask);
        
        _wf.set_to_buf(i%_param->nt_diag);
        //acc_vec[i] = _wf.acc();
        if(i%_param->nt_diag == 0){
            _wf.acc_buf();
            std::memcpy(&acc_vec[i],_wf.get_diag_buf(),_param->nt_diag*sizeof(cdouble));
            //norm_vec_idx = (int)(i/100);
            //norm_vec[norm_vec_idx] = _wf.norm();
            //std::cout<<"Ener: "<<(_ham.*(_ham.ener))(psi_row)<<"\n";
        }
    }

    std::cout<<"Norm: "<<_wf.norm()<<"\n";   
    std::string path = "results/end_psi2.dat";
    _wf.save_wf2(path);
    
    for(int i=0; i<_param->nt;i++){
        acc_vec[i] *= _accmask[i];
    }
    path = "results/acc.dat";
    write_array(acc_vec,_param->nt,path);
    path = "results/norm.dat";
    //write_array(norm_vec,norm_vec_size,path);
    delete acc_vec;
    //delete norm_vec;
    cdouble *temp;
    temp = new cdouble[_param->ni];
    for(int i=0; i<_param->ni; i++){
        temp[i] = _ham.get_dpotential()[i*_param->nk + 0];
    }
    path = "results/dpotential.dat";
    write_array(temp,_param->ni,path);
    for(int i=0; i<_param->ni; i++){
        temp[i] = _ham.get_potential()[i*_param->nk + 0];
    }
    path = "results/potential.dat";
    write_array(temp,_param->ni,path);
    delete temp;

    debug3("[TDSESolver->propagate] End propagate");
}

void TDSESolver::ipropagate_XZ(){
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
                psi_col[id_thread*nk + k] = _wf.get()[i*nk+k];
            (_ham.*(_ham.step_k))(&psi_col[id_thread*nk],0.0,0.0,i,1,id_thread);
            _wf.set_col(&psi_col[id_thread*nk],i);
        }
        
        #pragma omp parallel for schedule(dynamic)
        for(int k=0;k<nk;k++){
            int id_thread = omp_get_thread_num();
            for(int i=0;i<ni;i++)
                psi_row[id_thread*ni + i] = _wf.get()[i*nk+k];
            (_ham.*(_ham.step_i))(&psi_row[id_thread*ni],0.0,0.0,k,1,id_thread);
            _wf.set_row(&psi_row[id_thread*ni],k);
        }
        
        norm = _wf.norm();
        _wf /= norm;
        if(j%50==0){
            ener = (_ham.*(_ham.ener))(_wf.get());
            std::cout<<"Norm: "<< norm<<" Ener: "<<ener<<"\n";
        }
    }
    std::cout<<"Ener: "<<ener<<"\n";
    delete psi_col;
     delete psi_row;
}

void TDSESolver::propagate_XZ(){
    cdouble norm, ener;
    cdouble *acc_vec, *dip_vec;
    cdouble *psi_col, *psi_row;
    int idx;
    const int ni = _param->ni;
    const int nk = _param->nk;
    std::string path;

    acc_vec = new cdouble [_param->nt];
    dip_vec = new cdouble [_param->nt];
   
    psi_col = new cdouble [_param->nk*_param->n_threads];
    psi_row = new cdouble [_param->ni*_param->n_threads];

    //wf_ptr = _wf.get_buf();
    _wf.set_to_buf(0); 
    for(int j=0; j<_param->nt;j++){
        #pragma omp parallel for schedule(dynamic)
        for(int i=0;i<ni;i++){
            cdouble *wf_ptr;
            wf_ptr = _wf.get_buf();
            int id_thread = omp_get_thread_num();
            for(int k=0;k<nk;k++)
                psi_col[id_thread*nk + k] = wf_ptr[j%_param->nt_diag*ni*nk+i*nk+k];
            (_ham.*(_ham.step_k))(&psi_col[id_thread*nk],Afield_k[j],Bfield_k[j],i,0,id_thread);
            _wf.set_col_buf_mask(&psi_col[id_thread*nk],_kmask, i,(j+1)%_param->nt_diag);
        }
        
        #pragma omp parallel for schedule(dynamic)
        for(int k=0;k<nk;k++){
            cdouble *wf_ptr;
            wf_ptr = _wf.get_buf();
            int id_thread = omp_get_thread_num();
            for(int i=0;i<ni;i++)
                psi_row[id_thread*ni + i] = wf_ptr[(j+1)%_param->nt_diag*ni*nk+i*nk+k];
            (_ham.*(_ham.step_i))(&psi_row[id_thread*ni],Afield_i[j],Bfield_k[j],k,0,id_thread);
            _wf.set_row_buf_mask(&psi_row[id_thread*ni],_imask,k,(j+1)%_param->nt_diag);
        }
        //_wf.apply_mask_buf_RZ(_imask,_kmask,(j+1)%_param->nt_diag);
        //_wf.set_to_buf(j%_param->nt_diag);
        //acc_vec[j] = _wf.acc();
        //dip_vec[j] = _wf.dipole();
        if ((j+2)%_param->nt_diag==0 && j<(_param->nt-_param->nt%_param->nt_diag)){ 
            idx = j - _param->nt_diag + 2;
            _wf.dipole_buf();
            std::memcpy(&dip_vec[idx], _wf.get_diag_buf(), _param->nt_diag*sizeof(cdouble));
            _wf.acc_buf();
            std::memcpy(&acc_vec[idx], _wf.get_diag_buf(), _param->nt_diag*sizeof(cdouble));
 
            //norm = _wf.norm();
            //ener = (_ham.*(_ham.ener))(_wf.get());
            //std::cout<<j<<" ACC: "<< acc_vec[j]<<"\n";
        }
    }
    // Get last batch of diagnostics from last idx up to _paran.nt
    _wf.dipole_buf();
    std::memcpy(&dip_vec[idx+1],_wf.get_diag_buf(),(_param->nt%_param->nt_diag-1)*sizeof(cdouble));
    _wf.acc_buf();
    std::memcpy(&acc_vec[idx+1],_wf.get_diag_buf(),(_param->nt%_param->nt_diag-1)*sizeof(cdouble));

    cdouble valaccmask;
    for(int j=0; j<_param->nt; j++){
        valaccmask = _accmask[j];
        acc_vec[j] *= valaccmask;
        dip_vec[j] *= valaccmask;
    }
    write_array(acc_vec,_param->nt,_param->acc_path);
    write_array(dip_vec,_param->nt,_param->dip_path);

    delete acc_vec;
    delete dip_vec;
}


void TDSESolver::ipropagate_RZ(){
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
                psi_col[id_thread*nk + k] = _wf.get()[i*nk+k];
            (_ham.*(_ham.step_k))(&psi_col[id_thread*nk],0.0,0.0,i,1,id_thread);
            _wf.set_col(&psi_col[id_thread*nk],i);
        }
        
        #pragma omp parallel for schedule(dynamic)
        for(int k=0;k<nk;k++){
            int id_thread = omp_get_thread_num();
            for(int i=0;i<ni;i++)
                psi_row[id_thread*ni + i] = _wf.get()[i*nk+k];
            (_ham.*(_ham.step_i))(&psi_row[id_thread*ni],0.0,0.0,k,1,id_thread);
            _wf.set_row(&psi_row[id_thread*ni],k);
        }
        
        #pragma omp parallel for schedule(dynamic)
        for(int i=0;i<ni;i++){
            int id_thread = omp_get_thread_num();
            for(int k=0;k<nk;k++)
                psi_col[id_thread*nk + k] = _wf.get()[i*nk+k];
            (_ham.*(_ham.step_k))(&psi_col[id_thread*nk],0.0,0.0,i,1,id_thread);
            _wf.set_col(&psi_col[id_thread*nk],i);
        }

        norm = _wf.norm();
        _wf /= norm;
        if(j%50==0){
            ener = (_ham.*(_ham.ener))(_wf.get());
            std::cout<<"Norm: "<< norm<<" Ener: "<<ener<<"\n";
        }
    }
    std::cout<<"Ener: "<<ener<<"\n";
    delete psi_col;
     delete psi_row;
}

void TDSESolver::propagate_RZ(){
    cdouble norm, ener;
    cdouble *acc_vec, *dip_vec;
    cdouble *psi_col, *psi_row;
    int idx;
    const int ni = _param->ni;
    const int nk = _param->nk;
    std::string path;

    acc_vec = new cdouble [_param->nt];
    dip_vec = new cdouble [_param->nt];
   
    psi_col = new cdouble [_param->nk*_param->n_threads];
    psi_row = new cdouble [_param->ni*_param->n_threads];

    //wf_ptr = _wf.get_buf();
    _wf.set_to_buf(0); 
    for(int j=0; j<_param->nt;j++){
        #pragma omp parallel for schedule(dynamic)
        for(int i=0;i<ni;i++){
            cdouble *wf_ptr;
            wf_ptr = _wf.get_buf();
            int id_thread = omp_get_thread_num();
            for(int k=0;k<nk;k++)
                psi_col[id_thread*nk + k] = wf_ptr[j%_param->nt_diag*ni*nk+i*nk+k];
            (_ham.*(_ham.step_k))(&psi_col[id_thread*nk],Afield_k[j],Bfield_k[j],i,0,id_thread);
            _wf.set_col_buf_mask(&psi_col[id_thread*nk],_kmask, i,(j+1)%_param->nt_diag);
        }
        
        #pragma omp parallel for schedule(dynamic)
        for(int k=0;k<nk;k++){
            cdouble *wf_ptr;
            wf_ptr = _wf.get_buf();
            int id_thread = omp_get_thread_num();
            for(int i=0;i<ni;i++)
                psi_row[id_thread*ni + i] = wf_ptr[(j+1)%_param->nt_diag*ni*nk+i*nk+k];
            (_ham.*(_ham.step_i))(&psi_row[id_thread*ni],Afield_k[j],Bfield_k[j],k,0,id_thread);
            _wf.set_row_buf_mask(&psi_row[id_thread*ni],_imask,k,(j+1)%_param->nt_diag);
        }
        
        #pragma omp parallel for schedule(dynamic)
        for(int i=0;i<ni;i++){
            cdouble *wf_ptr;
            wf_ptr = _wf.get_buf();
            int id_thread = omp_get_thread_num();
            for(int k=0;k<nk;k++)
                psi_col[id_thread*nk + k] = wf_ptr[(j+1)%_param->nt_diag*ni*nk+i*nk+k];
            (_ham.*(_ham.step_k))(&psi_col[id_thread*nk],Afield_k[j],Bfield_k[j],i,0,id_thread);
            _wf.set_col_buf_mask(&psi_col[id_thread*nk],_kmask,i,(j+1)%_param->nt_diag);
        }
        //_wf.apply_mask_buf_RZ(_imask,_kmask,(j+1)%_param->nt_diag);
        //_wf.set_to_buf(j%_param->nt_diag);
        //acc_vec[j] = _wf.acc();
        //dip_vec[j] = _wf.dipole();
        if ((j+2)%_param->nt_diag==0 && j<(_param->nt-_param->nt%_param->nt_diag)){ 
            idx = j - _param->nt_diag + 2;
            _wf.dipole_buf();
            std::memcpy(&dip_vec[idx], _wf.get_diag_buf(), _param->nt_diag*sizeof(cdouble));
            _wf.acc_buf();
            std::memcpy(&acc_vec[idx], _wf.get_diag_buf(), _param->nt_diag*sizeof(cdouble));
 
            //norm = _wf.norm();
            //ener = (_ham.*(_ham.ener))(_wf.get());
            //std::cout<<j<<" Norm: "<< norm<<"\n";
        }
    }
    // Get last batch of diagnostics from last idx up to _paran.nt
    _wf.dipole_buf();
    std::memcpy(&dip_vec[idx+1],_wf.get_diag_buf(),(_param->nt%_param->nt_diag-1)*sizeof(cdouble));
    _wf.acc_buf();
    std::memcpy(&acc_vec[idx+1],_wf.get_diag_buf(),(_param->nt%_param->nt_diag-1)*sizeof(cdouble));

    cdouble valaccmask;
    for(int j=0; j<_param->nt; j++){
        valaccmask = _accmask[j];
        acc_vec[j] *= valaccmask;
        dip_vec[j] *= valaccmask;
    }
    //path = "results/acc.dat";
    write_array(acc_vec,_param->nt,_param->acc_path);
    //path = "results/dip.dat";
    write_array(dip_vec,_param->nt,_param->dip_path);

    delete acc_vec;
    delete dip_vec;
}


