#include <iostream>
#include "hamiltonian.h"

Hamiltonian::Hamiltonian(){
}

Hamiltonian::Hamiltonian(Parameters *param){
    _param = param;
    _ni = _param->ni;
    _nk = _param->nk;

    switch(_param->geometry){
        case X:
            _allocate_X();
            break;
        case XZ:
            _allocate_XZ();
            break;

        case RZ:
            _allocate_RZ();
            break;
    }    
}

void Hamiltonian::set_geometry(double *i, double *k, const double di, const double dk){
    _i = i; _k = k; _di = di; _dk = dk;
}

void Hamiltonian::set_potential(){

    switch(_param->geometry){
        case X:
            potential_X();
            break;
        case XZ:
            potential_XZ();
            break;
        case RZ:
            potential_RZ();
            break;
    }
        
}

void Hamiltonian::set_dpotential(){
    switch(_param->geometry){
        case X:
            dpotential_X();
            break;
        case XZ:
            dpotential_XZ();
            break;
        case RZ:
            dpotential_RZ();
            break;
    }
}

cdouble* Hamiltonian::get_potential(){
    return _potential;
}

cdouble* Hamiltonian::get_dpotential_i(){
    return _dpotential_i;
}

cdouble* Hamiltonian::get_dpotential_k(){
    return _dpotential_k;
}

void Hamiltonian::tridot(cdouble* aa, cdouble *bb, cdouble* cc, cdouble* vec, cdouble* out, const int n){
    // aa-> lower diagonal; bb-> main diagonal; cc-> upper diagonal
    out[0] = bb[0]*vec[0] + cc[0]*vec[1];
    out[n-1] = aa[n-1]*vec[n-2] + bb[n-1]*vec[n-1];
    for(int i=1;i<n-1;i++){
        out[i] = aa[i]*vec[i-1] + bb[i]*vec[i] + cc[i]*vec[i+1];
    }
}

void Hamiltonian::tdma(cdouble* aa, cdouble* bb, cdouble* cc, cdouble* dd, cdouble* out, const int n){
    // aa-> lower diagonal; bb-> main diagonal; cc-> upper diagonal
    
    cdouble wc;
    for(int i = 1; i<=n-1;i++){
        wc = aa[i]/bb[i-1];
        bb[i] = bb[i] - wc*cc[i-1];
        dd[i] = dd[i] - wc*dd[i-1];
    }

    out[n-1] = dd[n-1]/bb[n-1];
    
    for(int i = n-2;i>=0;i--){
        out[i] = (dd[i]-cc[i]*out[i+1])/bb[i];
    }
}



