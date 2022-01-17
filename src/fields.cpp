#include <iostream>
#include <cmath>
#include "fields.h"
#include "debug.h"
#include "parameters.h"

Field::Field(){
    _field = NULL;
}

Field::Field(double amp, double w, double phi, int env, double tmax,  double *t, int nt){
    _nt = nt; _tmax = tmax; _w = w; _t = t;
    _dt = _t[_nt-1]/(double)(_nt);
    _field = new double[_nt];

    for(int i=0; i<_nt; i++){
        if(env == 0){
            _field[i] = amp*env_sin2(_t[i])*sin(w*_t[i]+phi);
        }        
        else if(env == 1){
            _field[i] = amp*env_trap(_t[i])*sin(w*_t[i]+phi);
        }        
    }
    _flag = true;
}

void Field::calc_pot(){
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
    delete temp;
}

double Field::operator[](int i){
    return _field[i];
}

double Field::get(int i){
    return _field[i];
}

double* Field::get(){
    return _field;
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

