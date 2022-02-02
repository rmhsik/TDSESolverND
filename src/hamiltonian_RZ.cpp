#include <iostream>
#include "hamiltonian.h"

void Hamiltonian::_allocate_RZ(){
    step_i = &Hamiltonian::step_i_RZ;
    step_k = &Hamiltonian::step_k_RZ;
    ener = &Hamiltonian::ener_RZ;

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

}


void Hamiltonian::step_i_RZ(cdouble *psi_row, const int k, const int ti, const int imag, const int id_thread){
    cdouble Hr_du;
    cdouble Hr_d;
    cdouble Hr_dl;
    cdouble dt = imag==0 ? cdouble(1.0,0.0)*_param->dt: cdouble(0.0,-1.0)*_param->dt_ITP;

    cdouble bfield_k = (*Bfield_k)[ti];

    cdouble a = 1.0/(_di*_di);
    cdouble b = 1.0/(_di);
    cdouble c = 1.0/(2.0*_di); 
    for(int i=0;i<_ni;i++){
        Hr_du = -c*(b+0.5*1.0/(_i[i]));
        Hr_d  =  a + 0.5*_potential_fn(_i[i],_k[k],_t[ti]) + 0.5*0.125*bfield_k*bfield_k*_i[i]*_i[i];
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

void Hamiltonian::step_k_RZ(cdouble *psi_col, const int i, const int ti, const int imag, const int id_thread){
    cdouble Hz_du;
    cdouble Hz_d;
    cdouble Hz_dl;
    cdouble dt = imag==0 ? cdouble(1.0,0.0)*_param->dt: cdouble(0.0,-1.0)*_param->dt_ITP;
    
    cdouble afield_k = (*Afield_k)[ti];
    cdouble bfield_k = (*Bfield_k)[ti];

    cdouble a = 1.0/(_dk*_dk)+0.5*afield_k*afield_k/(C*C);
    cdouble b = 1.0/(2.0*_dk*_dk);
    Hz_du = -b - I*1.0/(2.0*C*_dk)*afield_k;
    Hz_dl = -b + I*1.0/(2.0*C*_dk)*afield_k;
    for(int k=0;k<_nk;k++){
        Hz_d  =  a + 0.5*_potential_fn(_i[i], _k[k],_t[ti])+ 0.5*0.125*bfield_k*bfield_k*_i[i]*_i[i];

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
            Hz_d[k]  =  a + 0.5*_potential_fn(_i[i],_k[k],0);
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
            Hr_d[i]  =  a + 0.5*_potential_fn(_i[i],_k[k],0);
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




