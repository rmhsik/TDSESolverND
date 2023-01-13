#include <iostream>
#include <cmath>
#include "fields.h"
#include "debug.h"
#include "parameters.h"
#include "utils.h"

Field::Field(){
    _field = NULL;
}

Field::Field(double amp, double w, double kw, double phi, int env, double tmax,  double *t, int nt){
    _amp = amp; _w = w; _kw = kw; _phi = phi;
    _nt = nt; _tmax = tmax; _w = w; _t = t;
    _dt = _t[_nt-1]/(double)(_nt);
    
    /*
    for(int i=0; i<_nt; i++){
        if(env == 0){
            _field[i] = amp*env_sin2(_t[i])*sin(w*_t[i]+phi);
        }        
        else if(env == 1){
            _field[i] = amp*env_trap(_t[i])*sin(w*_t[i]+phi);
        }        
    }
    _flag = true;
    */
}

void Field::set_params(Parameters *param){
    _param = param;

    _field = alloc2d<double>(_param->ni, _param->nk);
    _vecpot = alloc2d<double>(_param->ni, _param->nk);

    #pragma omp parallel for collapse(1) schedule(dynamic)
    for(int i=0; i<_param->ni; i++){
        for(int k=0; k<_param->nk; k++){
            _vecpot[i][k] = 0.0;
            _field[i][k]  = 0.0;
        }
    }
}

void Field::set_geometry(double *i, double *k, const double di, const double dk){
    _i = i; _k = k; _di = di; _dk = dk;
}

void Field::calc_field(double t){
    double env = env_sin2(t);
    #pragma omp parallel for collapse(1) schedule(dynamic) 
    for(int i = 0; i < _param->ni; i++){
        for(int k = 0; k < _param->nk; k++){
            _field[i][k] = env*plane_wave(_i[i], _k[k], t); 
        }
    }
}

void Field::calc_pot(){
    #pragma omp parallel for collapse(1) schedule(dynamic)
    for(int i=0;i<_param->ni; i++){
        for(int k=0; k<_param->nk; k++){
            _vecpot[i][k] += -C*_field[i][k]*_dt;
        }
    }

    /*
    double *temp;
    temp = new double[_nt];
    
    for(int i=0; i<_nt;i++){
        temp[i] = 0.0;
        for (int j=0; j<=i; j++){
            temp[i] += _field[j];
        }
        temp[i] *= -C*_dt;
    }
    for(int i=0; i<_nt;i++){
        _field[i] = temp[i];
    }
    delete[] temp;
    */
}

double Field::operator()(int i, int k){
    return _vecpot[i][k];
}

double Field::get(int i, int k){
    return _vecpot[i][k];
}

double** Field::get(){
    return _vecpot;
}

double Field::plane_wave(double i, double k, double t){
    double f = _amp*sin(_kw*i - _w*t);
    return f;
           
}

double Field::env_sin2(double ti){
    if (ti<_tmax){
        return pow(sin(M_PI*ti/_tmax),2);
    }
    else {
        return 0.0;
    }

}

double Field::env_trap(double ti){
    double T = 2*M_PI/_w;
    if (ti<_tmax+2.0*T){
         if (ti<T){
             return pow(sin(M_PI*ti/(2.0*T)),2);
         }
         else if (T<ti && ti<_tmax+T){
             return 1.0;
         }
         else{
             return pow(sin(M_PI*(ti-_tmax-T)/(2.0*T)+M_PI/2.0),2);
         }
    }
    else{
        return 0.0;
    }
}

Field::~Field(){
    delete[] _field;
    free2d<double>(&_field, _param->ni, _param->nk);
    free2d<double>(&_vecpot, _param->ni, _param->nk);
}
