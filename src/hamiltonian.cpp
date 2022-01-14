#include <iostream>
#include "hamiltonian.h"

Hamiltonian::Hamiltonian(){
}

Hamiltonian::Hamiltonian(Parameters param){
    _param = param;
    _ni = _param.ni;
     switch(_param.geometry){
            case X:
                step_i = &Hamiltonian::step_i_X;
                step_k = NULL;
                ener = &Hamiltonian::ener_X;
		break;
            case RZ:
                step_i = &Hamiltonian::step_i_RZ;
                step_k = &Hamiltonian::step_k_RZ;
                ener = &Hamiltonian::ener_RZ;
		break;
        }    
    _nk = _param.nk;
}

void Hamiltonian::set_geometry(double *i, double *k, const double di, const double dk){
    _i = i; _k = k; _di = di; _dk = dk;
}

void Hamiltonian::set_potential(){
    _potential = new cdouble[_ni*_nk];
    
    switch(_param.geometry){
        case X:
            potential_X();
            break;
        case RZ:
            potential_RZ();
            break;
    }
    
}

void Hamiltonian::set_dpotential(){
    _dpotential = new cdouble[_ni*_nk];
    switch(_param.geometry){
        case X:
            dpotential_X();
            break;
        case RZ:
            dpotential_RZ();
            break;
    }
}

void Hamiltonian::potential_X(){
    for(int i=0; i<_ni;i++){
        _potential[i*_nk + 0] = -1.0/sqrt(_i[i]*_i[i] + 2.0);
    }
}

void Hamiltonian::potential_RZ(){
    for(int i=0;i<_ni;i++){
        for(int k=0;k<_nk;k++){
            _potential[i*_nk + k] = -1.0/sqrt(_i[i]*_i[i] + _k[k]*_k[k]);
        }
    }
}

void Hamiltonian::dpotential_X(){
    _dpotential[0*_nk + 0] = (_potential[1*_nk + 0] - _potential[0*_nk + 0])/_di;
    _dpotential[(_ni-1)*_nk + 0] = (_potential[(_ni-1)*_nk + 0]-_potential[(_ni-2)*_nk + 0])/_di;
    for(int i=1; i<_ni-1;i++){
        _dpotential[i*_nk + 0] = (_potential[(i+1)*_nk + 0]-_potential[(i-1)*_nk + 0])/(2.0*_di);
    }
}

void Hamiltonian::dpotential_RZ(){
    //TODO
}

cdouble* Hamiltonian::get_potential(){
    return _potential;
}

cdouble* Hamiltonian::get_dpotential(){
    return _dpotential;
}

void Hamiltonian::step_i_X(cdouble *psi, double afield_i, double bfield_i, const int j, const int imag){
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
        H_du = -b + I*1.0/(2.0*C*_di)*afield_i;
        H_d  =  a + _potential[i*_nk + 0] + 0.5*afield_i*afield_i/(C*C);
        H_dl = -b - I*1.0/(2.0*C*_di)*afield_i;

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

cdouble Hamiltonian::ener_X(cdouble *psi){
    cdouble integral=0.0;

    cdouble *temp;
    cdouble *row;
    cdouble *H_du, *H_dl, *H_d; 
    temp = new cdouble[_ni];
    row = new cdouble[_ni];
    H_dl = new cdouble[_ni];
    H_d = new cdouble[_ni];
    H_du = new cdouble[_ni];
    
    cdouble a = 1.0/(_di*_di);
    cdouble b = 1.0/(2.0*_di*_di); 
    
    for(int i=0;i<_ni;i++){
        H_du[i] = -b;
        H_d[i]  =  a + _potential[i*_nk + 0];
        H_dl[i] = -b;
        row[i] = psi[i*_nk + 0];
    }
    tridot(H_dl, H_d, H_du, row, temp, _ni);
    
    for(int i=0; i<_ni; i++){
        integral += conj(row[i])*temp[i]*_di;
    }
    delete temp;
    delete row;
    delete H_du;
    delete H_d;
    delete H_dl; 
    return integral; 
}


void Hamiltonian::step_i_RZ(cdouble *psi_row, double afield_i, double bfield_i, const int k, const int imag ){
    cdouble Hr_du;
    cdouble Hr_d;
    cdouble Hr_dl;
    cdouble dt = imag==0 ? cdouble(1.0,0.0)*_param.dt: cdouble(0.0,-1.0)*_param.dt_ITP;

    cdouble *Mr_du, *Mpr_du;
    cdouble *Mr_d, *Mpr_d;
    cdouble *Mr_dl, *Mpr_dl;
    cdouble *lhs, *res;
     
    Mr_du = new cdouble[_ni];
    Mr_d  = new cdouble[_ni];
    Mr_dl = new cdouble[_ni];
    Mpr_du = new cdouble[_ni];
    Mpr_d  = new cdouble[_ni];
    Mpr_dl = new cdouble[_ni];
    lhs = new cdouble[_ni];
    res = new cdouble[_ni];

    cdouble a = 1.0/(_di*_di);
    cdouble b = 1.0/(_di);
    cdouble c = 1.0/(2.0*_di); 
    for(int i=0;i<_ni;i++){
        Hr_du = -c*(b+0.5*1.0/(_i[i]));
        Hr_d  =  a + 0.5*_potential[i*_nk + k];
        Hr_dl = -c*(b-0.5*1.0/(_i[i]));

        Mr_du[i] = I*Hr_du*dt/2.0;
        Mr_d[i]  = 1.0 + I*Hr_d*dt/2.0;
        Mr_dl[i] = I*Hr_dl*dt/2.0;
        Mpr_du[i] = -I*Hr_du*dt/2.0;
        Mpr_d[i]  = 1.0 - I*Hr_d*dt/2.0;
        Mpr_dl[i] = -I*Hr_dl*dt/2.0;
    }
    tridot(Mpr_du,Mpr_d,Mpr_dl,psi_row,lhs,_ni); 
    tdma(Mr_du,Mr_d,Mr_dl,lhs,res,_ni);
    for(int i=0; i<_ni; i++){
        psi_row[i] = res[i];
    }
    
    delete Mr_du;
    delete Mr_d;
    delete Mr_dl;
    delete Mpr_du;
    delete Mpr_d;
    delete Mpr_dl;
    delete lhs;
    delete res;
}

void Hamiltonian::step_k_RZ(cdouble *psi_col, double afield_k, double bfield_k, const int i, const int imag ){
    cdouble Hz_du;
    cdouble Hz_d;
    cdouble Hz_dl;
    cdouble dt = imag==0 ? cdouble(1.0,0.0)*_param.dt: cdouble(0.0,-1.0)*_param.dt_ITP;
    
    cdouble *Mz_du, *Mpz_du;
    cdouble *Mz_d, *Mpz_d;
    cdouble *Mz_dl, *Mpz_dl;
    cdouble *lhs, *res;
     
    Mz_du = new cdouble[_nk];
    Mz_d  = new cdouble[_nk];
    Mz_dl = new cdouble[_nk];
    Mpz_du = new cdouble[_nk];
    Mpz_d  = new cdouble[_nk];
    Mpz_dl = new cdouble[_nk];
    lhs = new cdouble[_nk];
    res = new cdouble[_nk];

    cdouble a = 1.0/(_dk*_dk);
    cdouble b = 1.0/(2.0*_dk*_dk);
    for(int k=0;k<_nk;k++){
        Hz_du = -b + I*1.0/(2.0*C*_dk)*afield_k;
        Hz_d  =  a + 0.5*_potential[i*_nk + k]+0.5*afield_k/(C*C);
        Hz_dl = -b - I*1.0/(2.0*C*_dk)*afield_k;
        Mz_du[k] = I*Hz_du*dt/4.0;
        Mz_d[k]  = 1.0 + I*Hz_d*dt/4.0;
        Mz_dl[k] = I*Hz_dl*dt/4.0;
        Mpz_du[k] = -I*Hz_du*dt/4.0;
        Mpz_d[k]  = 1.0 - I*Hz_d*dt/4.0;
        Mpz_dl[k] = -I*Hz_dl*dt/4.0;
    }
    tridot(Mpz_du,Mpz_d,Mpz_dl,psi_col,lhs,_nk); 
    tdma(Mz_du,Mz_d,Mz_dl,lhs,res,_nk);
    for(int k=0; k<_nk; k++){
        psi_col[k] = res[k];
    }
    
    delete Mz_du;
    delete Mz_d;
    delete Mz_dl;
    delete Mpz_du;
    delete Mpz_d;
    delete Mpz_dl;
    delete lhs;
    delete res;
}



cdouble Hamiltonian::ener_RZ(cdouble *psi){
    cdouble integral=0.0;

    cdouble *temp_r, *temp_z;
    cdouble *row, *col;
    cdouble *Hr_du, *Hr_dl, *Hr_d;
    cdouble *Hz_du, *Hz_dl, *Hz_d; 
    cdouble *wf_temp_r, *wf_temp_z;

    temp_r = new cdouble[_ni];
    temp_z = new cdouble[_nk];
    row = new cdouble[_ni];
    col = new cdouble[_nk];
    Hr_dl = new cdouble[_ni];
    Hr_d = new cdouble[_ni];
    Hr_du = new cdouble[_ni];
    Hz_dl = new cdouble[_nk];
    Hz_d = new cdouble[_nk];
    Hz_du = new cdouble[_nk];

    wf_temp_r = new cdouble[_ni*_nk];
    wf_temp_z = new cdouble[_ni*_nk];
    
    // Apply Hz to wavefunction
    {
    cdouble a = 1.0/(_dk*_dk);
    cdouble b = 1.0/(2.0*_dk*_dk); 
    for(int i=0;i<_ni;i++){
        for(int k=0;k<_nk;k++){
            Hz_du[k] = -b;
            Hz_d[k]  =  a + 0.5*_potential[i*_nk + k];
            Hz_dl[k] = -b;
            col[k] = psi[i*_nk + k];
        }
        tridot(Hz_dl, Hz_d, Hz_du, col, temp_z, _nk);

        for(int k=0; k<_nk; k++){
            wf_temp_z[i*_nk + k] = temp_z[k];
        }
    }
    }
    // Apply Hr to wavefunction
    {
    cdouble a = 1.0/(_di*_di);
    cdouble b = 1.0/(_di);
    cdouble c = 1.0/(2.0*_di);
    
    for(int k=0;k<_nk;k++){
        for(int i=0;i<_ni;i++){
            Hr_du[i] = -c*(b+0.5*1.0/(_i[i]));
            Hr_d[i]  =  a + 0.5*_potential[i*_nk + k];
            Hr_dl[i] = -c*(b-0.5*1.0/(_i[i]));;
            row[i] = psi[i*_nk + k];
        }
        tridot(Hr_dl, Hr_d, Hr_du, row, temp_r, _ni);

        for(int i=0; i<_ni; i++){
            wf_temp_r[i*_nk + k] = temp_r[i];
        }
    }
    }
    // Integrate 

    for(int i=0; i<_ni; i++){
        for(int k=0;k<_nk;k++)
            integral += 2.0*M_PI*_i[i]*conj(psi[i*_nk+k])*(wf_temp_z[i*_nk+k] + wf_temp_r[i*_nk+k])*_di*_dk;
    }
    delete temp_r;
    delete temp_z;
    delete row;
    delete col;
    delete Hr_du;
    delete Hr_d;
    delete Hr_dl;
    delete Hz_du;
    delete Hz_d;
    delete Hz_dl; 
    delete wf_temp_r;
    delete wf_temp_z;
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



