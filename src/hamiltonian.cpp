#include <iostream>
#include "hamiltonian.h"

Hamiltonian::Hamiltonian(){
}

Hamiltonian::Hamiltonian(Parameters param){
    _param = param;
    _ni = _param.ni;
    _nk = _param.nk;
}

void Hamiltonian::set_geometry(double *i, double *k, const double di, const double dk){
    _i = i; _k = k; _di = di; _dk = dk;
}

void Hamiltonian::set_potential(){
    _potential = new cdouble*[_ni];
    for (int i=0; i<_ni;i++){
        _potential[i] = new cdouble[_nk];
    }

    potential();
    
}

void Hamiltonian::set_dpotential(){
    _dpotential = new cdouble*[_ni];
    for(int i=0; i<_ni;i++){
        _dpotential[i] = new cdouble[_nk];
    }
    dpotential();
}

void Hamiltonian::potential(){
    for(int i=0; i<_ni;i++){
        for (int j=0; j<_nk;j++){
            _potential[i][j] = -1.0/sqrt(_i[i]*_i[i] + 2.0);
        }
    }
}

void Hamiltonian::dpotential(){
    _dpotential[0][0] = (_potential[1][0] - _potential[0][0])/_di;
    _dpotential[_ni-1][0] = (_potential[_ni-1][0]-_potential[_ni-2][0])/_di;
    for(int i=1; i<_ni-1;i++){
        _dpotential[i][0] = (_potential[i+1][0]-_potential[i-1][0])/(2.0*_di);
    }
}

cdouble** Hamiltonian::get_potential(){
    return _potential;
}

cdouble** Hamiltonian::get_dpotential(){
    return _dpotential;
}

void Hamiltonian::step_i(cdouble *psi, double afield_i, double bfield_i, const int imag){
    cdouble H_du;
    cdouble H_d;
    cdouble H_dl;
    cdouble dt = imag==0 ? cdouble(1.0,0.0)*_param.dt: cdouble(0.0,-1.0)*_param.dt_ITP;

    cdouble *M_du, *Mp_du;
    cdouble *M_d, *Mp_d;
    cdouble *M_dl, *Mp_dl;
    cdouble *lhs, *res;
     
    M_du = new cdouble[_ni];
    M_d  = new cdouble[_ni];
    M_dl = new cdouble[_ni];
    Mp_du = new cdouble[_ni];
    Mp_d  = new cdouble[_ni];
    Mp_dl = new cdouble[_ni];
    lhs = new cdouble[_ni];
    res = new cdouble[_ni];

    cdouble a = 1.0/(_di*_di);
    cdouble b = 1.0/(2.0*_di*_di); 
    
    for(int i=0;i<_ni;i++){
        H_du = -b + I*1.0/(2.0*c*_di)*afield_i;
        H_d  =  a + _potential[i][0] + 0.5*afield_i*afield_i/(c*c);
        H_dl = -b - I*1.0/(2.0*c*_di)*afield_i;

        M_du[i] = I*H_du*dt/2.0;
        M_d[i]  = 1.0 + I*H_d*dt/2.0;
        M_dl[i] = I*H_dl*dt/2.0;
        Mp_du[i] = -I*H_du*dt/2.0;
        Mp_d[i]  = 1.0 - I*H_d*dt/2.0;
        Mp_dl[i] = -I*H_dl*dt/2.0;
    }
    tridot(Mp_du,Mp_d,Mp_dl,psi,lhs,_ni); 
    tdma(M_du,M_d,M_dl,lhs,res,_ni);
    
    for(int i=0; i<_ni; i++){
        psi[i] = res[i];
    }
    
    delete M_du;
    delete M_d;
    delete M_dl;
    delete Mp_du;
    delete Mp_d;
    delete Mp_dl;
    delete lhs;
    delete res;
}

cdouble Hamiltonian::ener(cdouble *psi){
    cdouble integral=0.0;

    cdouble *temp;
    cdouble *H_du, *H_dl, *H_d; 
    temp = new cdouble[_ni];
    H_dl = new cdouble[_ni];
    H_d = new cdouble[_ni];
    H_du = new cdouble[_ni];
    
    cdouble a = 1.0/(_di*_di);
    cdouble b = 1.0/(2.0*_di*_di); 
    
    for(int i=0;i<_ni;i++){
        H_du[i] = -b;
        H_d[i]  =  a + _potential[i][0];
        H_dl[i] = -b;
    }
    tridot(H_dl, H_d, H_du, psi, temp, _ni);
    
    for(int i=0; i<_ni; i++){
        integral += conj(psi[i])*temp[i]*_di;
    }
    delete temp;
    delete H_du;
    delete H_d;
    delete H_dl; 
    return integral; 
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



