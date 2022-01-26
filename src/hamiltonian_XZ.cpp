#include <iostream>
#include "hamiltonian.h"

void Hamiltonian::_allocate_XZ(){
    step_i = &Hamiltonian::step_i_XZ;
    step_k = &Hamiltonian::step_k_XZ;
    ener = &Hamiltonian::ener_XZ;

    _Mk_du = new cdouble[_param->n_threads*_nk];
    _Mk_d  = new cdouble[_param->n_threads*_nk];
    _Mk_dl = new cdouble[_param->n_threads*_nk];
    _Mpk_du = new cdouble[_param->n_threads*_nk];
    _Mpk_d  = new cdouble[_param->n_threads*_nk];
    _Mpk_dl = new cdouble[_param->n_threads*_nk];
    _lhs_k = new cdouble[_param->n_threads*_nk];
    _res_k = new cdouble[_param->n_threads*_nk];
    
    _Mi_du = new cdouble[_param->n_threads*_ni];
    _Mi_d  = new cdouble[_param->n_threads*_ni];
    _Mi_dl = new cdouble[_param->n_threads*_ni];
    _Mpi_du = new cdouble[_param->n_threads*_ni];
    _Mpi_d  = new cdouble[_param->n_threads*_ni];
    _Mpi_dl = new cdouble[_param->n_threads*_ni];
    _lhs_i = new cdouble[_param->n_threads*_ni];
    _res_i = new cdouble[_param->n_threads*_ni];

    _potential = new cdouble[_ni*_nk];
    _dpotential_i = new cdouble[_ni*_nk];
    _dpotential_k = new cdouble[_ni*_nk];
}

void Hamiltonian::potential_XZ(){
    for(int i=0;i<_ni;i++){
        for(int k=0;k<_nk;k++){
            _potential[i*_nk + k] = -1.0/sqrt(_i[i]*_i[i] + _k[k]*_k[k] + 0.65);
        }
    }
}

void Hamiltonian::dpotential_XZ(){
    for(int i=0; i<_ni;i++){
        _dpotential_k[i*_nk + 0] = (_potential[i*_nk + 1] - _potential[i*_nk + 0])/_dk;
        _dpotential_k[i*_nk + _nk-1] = (_potential[i*_nk + _nk-1]-_potential[i*_nk + _nk-2])/_dk;
    }
    for(int i=0; i<_ni;i++){
        for(int k=1;k<_nk-1;k++)
            _dpotential_k[i*_nk + k] = (_potential[i*_nk + k+1]-_potential[i*_nk + k-1])/(2.0*_dk);
    }

    for(int k=0; k<_nk;k++){
        _dpotential_i[0*_nk + k] = (_potential[1*_nk + k] - _potential[0*_nk + k])/_di;
        _dpotential_i[(_ni-1)*_nk + k] = (_potential[(_ni-1)*_nk + k]-_potential[(_nk-2)*_nk + k])/_di;
    }
    for(int i=1; i<_ni-1;i++){
        for(int k=0;k<_nk;k++)
            _dpotential_i[i*_nk + k] = (_potential[(i+1)*_nk + k]-_potential[(i-1)*_nk + k])/(2.0*_di);
    }
}

void Hamiltonian::step_i_XZ(cdouble *psi_row, double afield_i, double bfield_k, const int k, const int imag, const int id_thread){
    cdouble Hx_du;
    cdouble Hx_d;
    cdouble Hx_dl;
    cdouble dt = imag==0 ? cdouble(1.0,0.0)*_param->dt: cdouble(0.0,-1.0)*_param->dt_ITP;


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
    cdouble dt = imag==0 ? cdouble(1.0,0.0)*_param->dt: cdouble(0.0,-1.0)*_param->dt_ITP;
    
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



