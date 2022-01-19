#include <iostream>
#include "hamiltonian.h"

Hamiltonian::Hamiltonian(){
}

Hamiltonian::Hamiltonian(Parameters param){
    _param = param;
    _ni = _param.ni;
    _nk = _param.nk;

    switch(_param.geometry){
        case X:
            step_i = &Hamiltonian::step_i_X;
            step_k = NULL;
            ener = &Hamiltonian::ener_X;
            _Mi_du = new cdouble[_ni];
            _Mi_d  = new cdouble[_ni];
            _Mi_dl = new cdouble[_ni];
            _Mpi_du = new cdouble[_ni];
            _Mpi_d  = new cdouble[_ni];
            _Mpi_dl = new cdouble[_ni];
            _lhs_i = new cdouble[_ni];
            _res_i = new cdouble[_ni];
            break;
        case XZ:
            step_i = &Hamiltonian::step_i_XZ;
            step_k = &Hamiltonian::step_k_XZ;
            ener = &Hamiltonian::ener_XZ;

            _Mk_du = new cdouble[n_threads*_nk];
            _Mk_d  = new cdouble[n_threads*_nk];
            _Mk_dl = new cdouble[n_threads*_nk];
            _Mpk_du = new cdouble[n_threads*_nk];
            _Mpk_d  = new cdouble[n_threads*_nk];
            _Mpk_dl = new cdouble[n_threads*_nk];
            _lhs_k = new cdouble[n_threads*_nk];
            _res_k = new cdouble[n_threads*_nk];
            
            _Mi_du = new cdouble[n_threads*_ni];
            _Mi_d  = new cdouble[n_threads*_ni];
            _Mi_dl = new cdouble[n_threads*_ni];
            _Mpi_du = new cdouble[n_threads*_ni];
            _Mpi_d  = new cdouble[n_threads*_ni];
            _Mpi_dl = new cdouble[n_threads*_ni];
            _lhs_i = new cdouble[n_threads*_ni];
            _res_i = new cdouble[n_threads*_ni];
            break;

        case RZ:
            step_i = &Hamiltonian::step_i_RZ;
            step_k = &Hamiltonian::step_k_RZ;
            ener = &Hamiltonian::ener_RZ;

            _Mk_du = new cdouble[n_threads*_nk];
            _Mk_d  = new cdouble[n_threads*_nk];
            _Mk_dl = new cdouble[n_threads*_nk];
            _Mpk_du = new cdouble[n_threads*_nk];
            _Mpk_d  = new cdouble[n_threads*_nk];
            _Mpk_dl = new cdouble[n_threads*_nk];
            _lhs_k = new cdouble[n_threads*_nk];
            _res_k = new cdouble[n_threads*_nk];
            
            _Mi_du = new cdouble[n_threads*_ni];
            _Mi_d  = new cdouble[n_threads*_ni];
            _Mi_dl = new cdouble[n_threads*_ni];
            _Mpi_du = new cdouble[n_threads*_ni];
            _Mpi_d  = new cdouble[n_threads*_ni];
            _Mpi_dl = new cdouble[n_threads*_ni];
            _lhs_i = new cdouble[n_threads*_ni];
            _res_i = new cdouble[n_threads*_ni];
            break;
        }    
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
        case XZ:
            potential_XZ();
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
        case XZ:
            dpotential_XZ();
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

void Hamiltonian::potential_XZ(){
    for(int i=0;i<_ni;i++){
        for(int k=0;k<_nk;k++){
            _potential[i*_nk + k] = -1.0/sqrt(_i[i]*_i[i] + _k[k]*_k[k] + 0.65);
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
    for(int i=0; i<_ni;i++){
        _dpotential[i*_nk + 0] = (_potential[i*_nk + 1] - _potential[i*_nk + 0])/_dk;
        _dpotential[i*_nk + _nk-1] = (_potential[i*_nk + _nk-1]-_potential[i*_nk + _nk-2])/_dk;
    }
    for(int i=0; i<_ni;i++){
        for(int k=0;k<_nk-1;k++)
            _dpotential[i*_nk + k] = (_potential[i*_nk + k+1]-_potential[i*_nk + k-1])/(2.0*_dk);
    }
}

void Hamiltonian::dpotential_XZ(){
    for(int i=0; i<_ni;i++){
        _dpotential[i*_nk + 0] = (_potential[i*_nk + 1] - _potential[i*_nk + 0])/_dk;
        _dpotential[i*_nk + _nk-1] = (_potential[i*_nk + _nk-1]-_potential[i*_nk + _nk-2])/_dk;
    }
    for(int i=0; i<_ni;i++){
        for(int k=0;k<_nk-1;k++)
            _dpotential[i*_nk + k] = (_potential[i*_nk + k+1]-_potential[i*_nk + k-1])/(2.0*_dk);
    }
}

cdouble* Hamiltonian::get_potential(){
    return _potential;
}

cdouble* Hamiltonian::get_dpotential(){
    return _dpotential;
}

void Hamiltonian::step_i_X(cdouble *psi, double afield_i, double bfield_i, const int j, const int imag,const int id_thread){
    cdouble H_du;
    cdouble H_d;
    cdouble H_dl;
    cdouble dt = imag==0 ? cdouble(1.0,0.0)*_param.dt: cdouble(0.0,-1.0)*_param.dt_ITP;

    cdouble a = 1.0/(_di*_di);
    cdouble b = 1.0/(2.0*_di*_di); 
    for(int i=0;i<_ni;i++){
        H_du = -b + I*1.0/(2.0*C*_di)*afield_i;
        H_d  =  a + _potential[i*_nk + 0] + 0.5*afield_i*afield_i/(C*C);
        H_dl = -b - I*1.0/(2.0*C*_di)*afield_i;

        _Mi_du[i] = I*H_du*dt/2.0;
        _Mi_d[i]  = 1.0 + I*H_d*dt/2.0;
        _Mi_dl[i] = I*H_dl*dt/2.0;
        _Mpi_du[i] = -I*H_du*dt/2.0;
        _Mpi_d[i]  = 1.0 - I*H_d*dt/2.0;
        _Mpi_dl[i] = -I*H_dl*dt/2.0;
    }
    tridot(_Mpi_du,_Mpi_d,_Mpi_dl,psi,_lhs_i,_ni); 
    tdma(_Mi_du,_Mi_d,_Mi_dl,_lhs_i,_res_i,_ni);
    for(int i=0; i<_ni; i++){
        psi[i] = _res_i[i];
    }
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

void Hamiltonian::step_i_XZ(cdouble *psi_row, double afield_i, double bfield_k, const int k, const int imag, const int id_thread){
    cdouble Hx_du;
    cdouble Hx_d;
    cdouble Hx_dl;
    cdouble dt = imag==0 ? cdouble(1.0,0.0)*_param.dt: cdouble(0.0,-1.0)*_param.dt_ITP;


    cdouble a = 1.0/(_di*_di) + 0.5*afield_i*afield_i/(C*C);
    cdouble b = 1.0/(2.0*_di*_di);
    Hx_du = -b - I*1.0/(2.0*C*_di)*afield_i;
    Hx_dl = -b + I*1.0/(2.0*C*_di)*afield_i;
    for(int i=0;i<_ni;i++){
        Hx_d  =  a + 0.5*_potential[i*_nk + k] + 0.5*0.125*bfield_k*bfield_k*_i[i]*_i[i];

        _Mi_du[id_thread*_ni + i] = I*Hx_du*dt/2.0;
        _Mi_d[id_thread*_ni + i]  = 1.0 + I*Hx_d*dt/2.0;
        _Mi_dl[id_thread*_ni + i] = I*Hx_dl*dt/2.0;
        _Mpi_du[id_thread*_ni + i] = -I*Hx_du*dt/2.0;
        _Mpi_d[id_thread*_ni + i]  = 1.0 - I*Hx_d*dt/2.0;
        _Mpi_dl[id_thread*_ni + i] = -I*Hx_dl*dt/2.0;
    }
    tridot(&_Mpi_dl[id_thread*_ni],
           &_Mpi_d[id_thread*_ni],
           &_Mpi_du[id_thread*_ni],
           psi_row,&_lhs_i[id_thread*_ni],_ni); 
    tdma(&_Mi_dl[id_thread*_ni],
         &_Mi_d[id_thread*_ni],
         &_Mi_du[id_thread*_ni],
         &_lhs_i[id_thread*_ni],&_res_i[id_thread*_ni],_ni);
    for(int i=0; i<_ni; i++){
        psi_row[i] = _res_i[id_thread*_ni + i];
    }
}

void Hamiltonian::step_k_XZ(cdouble *psi_col, double afield_k, double bfield_k, const int i, const int imag, const int id_thread){
    cdouble Hz_du;
    cdouble Hz_d;
    cdouble Hz_dl;
    cdouble dt = imag==0 ? cdouble(1.0,0.0)*_param.dt: cdouble(0.0,-1.0)*_param.dt_ITP;
    
    cdouble a = 1.0/(_dk*_dk)+0.5*afield_k*afield_k/(C*C);
    cdouble b = 1.0/(2.0*_dk*_dk);
    Hz_du = -b - I*1.0/(2.0*C*_dk)*afield_k;
    Hz_dl = -b + I*1.0/(2.0*C*_dk)*afield_k;
    for(int k=0;k<_nk;k++){
        Hz_d  =  a + 0.5*_potential[i*_nk + k]+ 0.5*0.125*bfield_k*bfield_k*_i[i]*_i[i];

        _Mk_du[id_thread*_nk + k] = I*Hz_du*dt/2.0;
        _Mk_d[id_thread*_nk + k]  = 1.0 + I*Hz_d*dt/2.0;
        _Mk_dl[id_thread*_nk + k] = I*Hz_dl*dt/2.0;
        _Mpk_du[id_thread*_nk + k] = -I*Hz_du*dt/2.0;
        _Mpk_d[id_thread*_nk + k]  = 1.0 - I*Hz_d*dt/2.0;
        _Mpk_dl[id_thread*_nk + k] = -I*Hz_dl*dt/2.0;
    }
    tridot(&_Mpk_dl[id_thread*_nk],
           &_Mpk_d[id_thread*_nk],
           &_Mpk_du[id_thread*_nk],
           psi_col, &_lhs_k[id_thread*_nk], _nk); 
    tdma(&_Mk_dl[id_thread*_nk],
         &_Mk_d[id_thread*_nk],
         &_Mk_du[id_thread*_nk],
         &_lhs_k[id_thread*_nk],&_res_k[id_thread*_nk],_nk);
    for(int k=0; k<_nk; k++){
        psi_col[k] = _res_k[id_thread*_nk + k];
    }
}

cdouble Hamiltonian::ener_XZ(cdouble *psi){
    cdouble integral=0.0;

    cdouble *temp_x, *temp_z;
    cdouble *row, *col;
    cdouble *Hx_du, *Hx_dl, *Hx_d;
    cdouble *Hz_du, *Hz_dl, *Hz_d; 
    cdouble *wf_temp_x, *wf_temp_z;

    temp_x = new cdouble[_ni];
    temp_z = new cdouble[_nk];
    row = new cdouble[_ni];
    col = new cdouble[_nk];
    Hx_dl = new cdouble[_ni];
    Hx_d = new cdouble[_ni];
    Hx_du = new cdouble[_ni];
    Hz_dl = new cdouble[_nk];
    Hz_d = new cdouble[_nk];
    Hz_du = new cdouble[_nk];

    wf_temp_x = new cdouble[_ni*_nk];
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
    // Apply Hx to wavefunction
    {
    cdouble a = 1.0/(_di*_di);
    cdouble b = 1.0/(2.0*_di*_di); 
    for(int k=0;k<_nk;k++){
        for(int i=0;i<_ni;i++){
            Hx_du[i] = -b;
            Hx_d[i]  =  a + 0.5*_potential[i*_nk + k];
            Hx_dl[i] = -b;
            row[i] = psi[i*_nk + k];
        }
        tridot(Hx_dl, Hx_d, Hx_du, row, temp_x, _ni);

        for(int i=0; i<_ni; i++){
            wf_temp_x[i*_nk + k] = temp_x[i];
        }
    }
    }
    // Integrate 

    for(int i=0; i<_ni; i++){
        for(int k=0;k<_nk;k++)
            integral += conj(psi[i*_nk+k])*(wf_temp_z[i*_nk+k] + wf_temp_x[i*_nk+k])*_di*_dk;
    }
    delete temp_x;
    delete temp_z;
    delete row;
    delete col;
    delete Hx_du;
    delete Hx_d;
    delete Hx_dl;
    delete Hz_du;
    delete Hz_d;
    delete Hz_dl; 
    delete wf_temp_x;
    delete wf_temp_z;
    return integral; 
}


void Hamiltonian::step_i_RZ(cdouble *psi_row, double afield_i, double bfield_k, const int k, const int imag, const int id_thread){
    cdouble Hr_du;
    cdouble Hr_d;
    cdouble Hr_dl;
    cdouble dt = imag==0 ? cdouble(1.0,0.0)*_param.dt: cdouble(0.0,-1.0)*_param.dt_ITP;


    cdouble a = 1.0/(_di*_di);
    cdouble b = 1.0/(_di);
    cdouble c = 1.0/(2.0*_di); 
    for(int i=0;i<_ni;i++){
        Hr_du = -c*(b+0.5*1.0/(_i[i]));
        Hr_d  =  a + 0.5*_potential[i*_nk + k] + 0.5*0.125*bfield_k*bfield_k*_i[i]*_i[i];
        Hr_dl = -c*(b-0.5*1.0/(_i[i]));

        _Mi_du[id_thread*_ni + i] = I*Hr_du*dt/2.0;
        _Mi_d[id_thread*_ni + i]  = 1.0 + I*Hr_d*dt/2.0;
        _Mi_dl[id_thread*_ni + i] = I*Hr_dl*dt/2.0;
        _Mpi_du[id_thread*_ni + i] = -I*Hr_du*dt/2.0;
        _Mpi_d[id_thread*_ni + i]  = 1.0 - I*Hr_d*dt/2.0;
        _Mpi_dl[id_thread*_ni + i] = -I*Hr_dl*dt/2.0;
    }
    tridot(&_Mpi_dl[id_thread*_ni],
           &_Mpi_d[id_thread*_ni],
           &_Mpi_du[id_thread*_ni],
           psi_row,&_lhs_i[id_thread*_ni],_ni); 
    tdma(&_Mi_dl[id_thread*_ni],
         &_Mi_d[id_thread*_ni],
         &_Mi_du[id_thread*_ni],
         &_lhs_i[id_thread*_ni],&_res_i[id_thread*_ni],_ni);
    for(int i=0; i<_ni; i++){
        psi_row[i] = _res_i[id_thread*_ni + i];
    }
}

void Hamiltonian::step_k_RZ(cdouble *psi_col, double afield_k, double bfield_k, const int i, const int imag, const int id_thread){
    cdouble Hz_du;
    cdouble Hz_d;
    cdouble Hz_dl;
    cdouble dt = imag==0 ? cdouble(1.0,0.0)*_param.dt: cdouble(0.0,-1.0)*_param.dt_ITP;
    
    cdouble a = 1.0/(_dk*_dk)+0.5*afield_k*afield_k/(C*C);
    cdouble b = 1.0/(2.0*_dk*_dk);
    Hz_du = -b - I*1.0/(2.0*C*_dk)*afield_k;
    Hz_dl = -b + I*1.0/(2.0*C*_dk)*afield_k;
    for(int k=0;k<_nk;k++){
        Hz_d  =  a + 0.5*_potential[i*_nk + k]+ 0.5*0.125*bfield_k*bfield_k*_i[i]*_i[i];

        _Mk_du[id_thread*_nk + k] = I*Hz_du*dt/4.0;
        _Mk_d[id_thread*_nk + k]  = 1.0 + I*Hz_d*dt/4.0;
        _Mk_dl[id_thread*_nk + k] = I*Hz_dl*dt/4.0;
        _Mpk_du[id_thread*_nk + k] = -I*Hz_du*dt/4.0;
        _Mpk_d[id_thread*_nk + k]  = 1.0 - I*Hz_d*dt/4.0;
        _Mpk_dl[id_thread*_nk + k] = -I*Hz_dl*dt/4.0;
    }
    tridot(&_Mpk_dl[id_thread*_nk],
           &_Mpk_d[id_thread*_nk],
           &_Mpk_du[id_thread*_nk],
           psi_col, &_lhs_k[id_thread*_nk], _nk); 
    tdma(&_Mk_dl[id_thread*_nk],
         &_Mk_d[id_thread*_nk],
         &_Mk_du[id_thread*_nk],
         &_lhs_k[id_thread*_nk],&_res_k[id_thread*_nk],_nk);
    for(int k=0; k<_nk; k++){
        psi_col[k] = _res_k[id_thread*_nk + k];
    }
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



