#include <iostream>
#include <tuple>
#include <string>
#include "tdsesolver.h"
#include "debug.h"
#include "utils.h"

TDSESolver::TDSESolver(){
}

TDSESolver::TDSESolver(Parameters param){
    std::string path;
    _param = param;
    setup_time();

    setup_geometry();
    setup_fields();
    setup_wf();
    setup_masks();
    }

void TDSESolver::setup_time(){
    std::string path;
    std::tie(_t,_dt) = linspace<double>(0.0,_param.tmax_sim, _param.nt);
    path = "results/time.dat";
    write_array(_t,_param.nt,path);
}

void TDSESolver::setup_geometry(){
    std::string path;
    switch(_param.geometry){
        case X:
            std::tie(_i,_di) = linspace<double>(_param.imin,_param.imax,_param.ni);
            std::tie(_k,_dk) = linspace<double>(_param.kmin,_param.kmax,_param.nk);
            break;
        case XZ:
            std::tie(_i,_di) = linspace<double>(_param.imin,_param.imax,_param.ni);
            std::tie(_k,_dk) = linspace<double>(_param.kmin,_param.kmax,_param.nk);
            break;
        case RZ:
            std::tie(_i,_di) = linspace<double>(0.0,_param.imax,_param.ni);
            std::tie(_k,_dk) = linspace<double>(_param.kmin,_param.kmax, _param.nk);
            
            for(int i=0; i<_param.ni;i++){
                _i[i] += 0.5;
            }
            break;
    }
    path = "results/i.dat";
    write_array(_i,_param.ni,path);
    path = "results/k.dat";
    write_array(_k,_param.nk,path);
}

void TDSESolver::setup_fields(){
    std::string path;
    switch(_param.geometry){
        case X:
            Afield_i = Field(_param.E0i, _param.w0Ei, _param.phiEi, _param.env, _param.tmax_ev,_t, _param.nt);
            Afield_i.calc_pot();
            path = "results/Afield_i.dat";
            write_array(Afield_i.get(),_param.nt,path);
            break;

        case XZ:
            Afield_i = Field(_param.E0i, _param.w0Ei, _param.phiEi, _param.env, _param.tmax_ev, _t, _param.nt);
            Afield_k = Field(_param.E0k, _param.w0Ek, _param.phiEk, _param.env, _param.tmax_ev, _t, _param.nt);
            Bfield_i = Field(_param.B0i, _param.w0Bi, _param.phiBi, _param.env, _param.tmax_ev, _t, _param.nt);
            Bfield_k = Field(_param.B0k, _param.w0Bk, _param.phiBk, _param.env, _param.tmax_ev, _t, _param.nt);
            Afield_i.calc_pot();
            Afield_k.calc_pot();
            path = "results/Afield_i.dat";
            write_array(Afield_i.get(),_param.nt,path);
            path = "results/Afield_k.dat";
            write_array(Afield_k.get(),_param.nt,path);
            path = "results/Bfield_i.dat";
            write_array(Bfield_i.get(),_param.nt,path);
            path = "results/Bfield_k.dat";
            write_array(Bfield_k.get(),_param.nt,path);
            break;

        case RZ:
            Afield_k = Field(_param.E0k, _param.w0Ek, _param.phiEk, _param.env, _param.tmax_ev, _t, _param.nt);
            Bfield_k = Field(_param.B0k, _param.w0Bk, _param.phiBk, _param.env, _param.tmax_ev, _t, _param.nt);
            Afield_k.calc_pot();
            path = "results/Afield_k.dat";
            write_array(Afield_k.get(),_param.nt,path);
            path = "results/Bfield_k.dat";
            write_array(Bfield_k.get(),_param.nt,path);
            break;
    }
}

void TDSESolver::setup_wf(){
    _wf = WF(_param, _i, _k, _di, _dk);
    _wf.gaussian(0.0,0.0,1.0);
    cdouble norm = _wf.norm();
    std::cout<<norm<<std::endl;
    _wf /= norm;
    norm = _wf.norm();
    std::cout<<norm<<std::endl;

    std::string path = "results/init_psi2.dat";
    _wf.save_wf2(path);
}

void TDSESolver::setup_masks(){
    _accmask = new cdouble[_param.nt];
    _imask = new cdouble[_param.ni];
    _kmask = new cdouble[_param.nk];
    double tmax_mask = 0.0;
    double ib = _param.imax/10.0;
    double kb = _param.kmax/10.0;
    double gamma = 1.0;
    std::string path;
    
    if (_param.env==SIN2){
        tmax_mask = _param.tmax_ev;
    }
    else{
        tmax_mask = _param.tmax_ev + _param.period*2.0;
    }

    for(int i=0; i<_param.nt;i++){
        if(_t[i]<tmax_mask){
            _accmask[i] = 1.0;
        }
        else{
            _accmask[i] = exp(-pow(_t[i]-tmax_mask,2)/0.5);
        }
    }
    path = "results/accmask.dat";
    write_array(_accmask,_param.nt,path);

    switch(_param.geometry){
        case X:
           for(int i=0; i<_param.ni;i++){
                if(_i[i]<_i[0]+ib){
                    _imask[i] = pow(cos(M_PI*(_i[i] - (_i[0]+ib))*gamma/(2*ib)),1.0/8.0);
                }
                else if (_i[i]>(_i[_param.ni-1]-ib)){
                    _imask[i] = pow(cos(M_PI*(_i[i] - (_i[_param.ni-1]-ib))*gamma/(2*ib)),1.0/8.0);
                }
                else{
                    _imask[i] = 1.0;
                }
            }
            _kmask[0] = 1.0;
            path = "results/imask.dat";
            write_array(_imask,_param.ni,path);
            path = "results/kmask.dat";
            write_array(_kmask,_param.nk,path);

            break;

        case XZ:
            for(int i=0; i<_param.ni;i++){
                if(_i[i]<_i[0]+ib){
                    _imask[i] = pow(cos(M_PI*(_i[i] - (_i[0]+ib))*gamma/(2*ib)),1.0/8.0);
                }
                else if (_i[i]>(_i[_param.ni-1]-ib)){
                    _imask[i] = pow(cos(M_PI*(_i[i] - (_i[_param.ni-1]-ib))*gamma/(2*ib)),1.0/8.0);
                }
                else{
                    _imask[i] = 1.0;
                }
            }

            for(int i=0; i<_param.nk;i++){
                if(_k[i]<_k[0]+kb){
                    _kmask[i] = pow(cos(M_PI*(_k[i] - (_k[0]+kb))*gamma/(2*kb)),1.0/8.0);
                }
                else if (_k[i]>(_k[_param.nk-1]-kb)){
                    _kmask[i] = pow(cos(M_PI*(_k[i] - (_k[_param.nk-1]-kb))*gamma/(2*kb)),1.0/8.0);
                }
                else{
                    _kmask[i] = 1.0;
                }
            }
            path = "results/imask.dat";
            write_array(_imask,_param.ni,path);
            path = "results/kmask.dat";
            write_array(_kmask,_param.nk,path);

            break; 

        case RZ:
            for(int i=0; i<_param.ni;i++){
                if (_i[i]>(_i[_param.ni-1]-ib)){
                    _imask[i] = pow(cos(M_PI*(_i[i] - (_i[_param.ni-1]-ib))*gamma/(2*ib)),1.0/8.0);
                }
                else{
                    _imask[i] = 1.0;
                }
            }

            for(int i=0; i<_param.nk;i++){
                if(_k[i]<_k[0]+kb){
                    _kmask[i] = pow(cos(M_PI*(_k[i] - (_k[0]+kb))*gamma/(2*kb)),1.0/8.0);
                }
                else if (_k[i]>(_k[_param.nk-1]-kb)){
                    _kmask[i] = pow(cos(M_PI*(_k[i] - (_k[_param.nk-1]-kb))*gamma/(2*kb)),1.0/8.0);
                }
                else{
                    _kmask[i] = 1.0;
                }
            }
            path = "results/imask.dat";
            write_array(_imask,_param.ni,path);
            path = "results/kmask.dat";
            write_array(_kmask,_param.nk,path);

            break;
    }

}
