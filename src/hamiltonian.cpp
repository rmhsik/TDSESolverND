#include <iostream>
#include "hamiltonian.h"
#include "potential.h"

Hamiltonian::Hamiltonian(){
}

Hamiltonian::Hamiltonian(Parameters *param){
    _param = param;
    _ni = _param->ni;
    _nk = _param->nk;
    _nt = _param->nt;

    switch(_param->geometry){
        case X:
            _allocate_X();
            if(_param->use_potential == 0)
                _potential = &potential_X;
            break;
        case XZ:
            _allocate_XZ();
            if(_param->use_potential == 0)
                _potential = &potential_XZ;
            break;

        case RZ:
            _allocate_RZ();
            if(_param->use_potential == 0)
                _potential = &potential_RZ;
            break;
    }    
    if(_param->use_potential != 0 )
            _potential = &potential;
}

void Hamiltonian::set_geometry(double *i, double *k, double *t, const double di, const double dk, const double dt){
    _i = i; _k = k; _t = t; _di = di; _dk = dk; _dt = dt;
}

void Hamiltonian::set_fields(Field* field1, Field* field2, Field* field3, Field* field4){
    Afield_i = field1;
    Afield_k = field2;
    Bfield_i = field3;
    Bfield_k = field4;

}

cdouble Hamiltonian::_potential_fn(double i, double k, double t){
    return (*_potential)(i,k,t,this);
    //return potential(i,k,t, this);
}

cdouble Hamiltonian::dpotential_i(double i, double k){
    return (_potential_fn(i + _di,k,0) - _potential_fn(i - _di,k,0))/(2.0*_di);
}

cdouble Hamiltonian::dpotential_k(double i, double k){
    return (_potential_fn(i,k + _dk,0) - _potential_fn(i,k - _dk,0))/(2.0*_dk);
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

Hamiltonian::~Hamiltonian(){
    delete[] _Mk_du;
    delete[] _Mk_d;
    delete[] _Mk_dl;
    delete[] _Mpk_du;
    delete[] _Mpk_d;
    delete[] _Mpk_dl;
    delete[] _Mi_du;
    delete[] _Mi_d;
    delete[] _Mi_dl;
    delete[] _Mpi_du;
    delete[] _Mpi_d;
    delete[] _Mpi_dl;
    delete[] _lhs_i;
    delete[] _lhs_k;
    delete[] _res_i;
    delete[] _res_k;
}
