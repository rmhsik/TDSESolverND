#include <iostream>
#include "hamiltonian.h"
#include <omp.h>
void Hamiltonian::_allocate_XYZ(){
    step_i = &Hamiltonian::step_i_XYZ;
    step_k = &Hamiltonian::step_k_XYZ;
    step_j = &Hamiltonian::step_j_XYZ;
    ener = &Hamiltonian::ener_XYZ;

    _Mk_du = new cdouble[_param->n_threads*_nk];
    _Mk_d  = new cdouble[_param->n_threads*_nk];
    _Mk_dl = new cdouble[_param->n_threads*_nk];
    _Mpk_du = new cdouble[_param->n_threads*_nk];
    _Mpk_d  = new cdouble[_param->n_threads*_nk];
    _Mpk_dl = new cdouble[_param->n_threads*_nk];
    _lhs_k = new cdouble[_param->n_threads*_nk];
    _res_k = new cdouble[_param->n_threads*_nk];

    _Mj_du = new cdouble[_param->n_threads*_nj];
    _Mj_d  = new cdouble[_param->n_threads*_nj];
    _Mj_dl = new cdouble[_param->n_threads*_nj];
    _Mpj_du = new cdouble[_param->n_threads*_nj];
    _Mpj_d  = new cdouble[_param->n_threads*_nj];
    _Mpj_dl = new cdouble[_param->n_threads*_nj];
    _lhs_j = new cdouble[_param->n_threads*_nj];
    _res_j = new cdouble[_param->n_threads*_nj];

    _Mi_du = new cdouble[_param->n_threads*_ni];
    _Mi_d  = new cdouble[_param->n_threads*_ni];
    _Mi_dl = new cdouble[_param->n_threads*_ni];
    _Mpi_du = new cdouble[_param->n_threads*_ni];
    _Mpi_d  = new cdouble[_param->n_threads*_ni];
    _Mpi_dl = new cdouble[_param->n_threads*_ni];
    _lhs_i = new cdouble[_param->n_threads*_ni];
    _res_i = new cdouble[_param->n_threads*_ni];
}

void Hamiltonian::step_i_XYZ(cdouble *psi_i_row, const int j, const int k, const int ti, const int imag, const int id_thread){
    cdouble Hx_du;
    cdouble Hx_d;
    cdouble Hx_dl;
    cdouble dt = imag==0 ? cdouble(1.0,0.0)*_param->dt: cdouble(0.0,-1.0)*_param->dt_ITP;

    cdouble afield_i = (*Afield_i)[ti];
    cdouble bfield_k = (*Bfield_k)[ti];

    cdouble a = 1.0/(_di*_di) + 0.5*afield_i*afield_i/(C*C);
    cdouble b = 1.0/(2.0*_di*_di);
    Hx_du = -b - I*1.0/(2.0*C*_di)*afield_i + I/(2.0*_di)*bfield_k*_j[j];
    Hx_dl = -b + I*1.0/(2.0*C*_di)*afield_i - I/(2.0*_di)*bfield_k*_j[j];
    for(int i=0;i<_ni;i++){
        Hx_d  =  a + 1.0/3.0*_potential_fn(_i[i],_j[j],_k[k],_t[ti]) + 1.0/3.0*0.125*bfield_k*bfield_k*(_i[i]*_i[i]+_j[j]*_j[j]);

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
           psi_i_row,&_lhs_i[id_thread*_ni],_ni); 
    tdma(&_Mi_dl[id_thread*_ni],
         &_Mi_d[id_thread*_ni],
         &_Mi_du[id_thread*_ni],
         &_lhs_i[id_thread*_ni],&_res_i[id_thread*_ni],_ni);
    for(int i=0; i<_ni; i++){
        psi_i_row[i] = _res_i[id_thread*_ni + i];
    }
}

void Hamiltonian::step_j_XYZ(cdouble *psi_j_row, const int i, const int k, const int ti, const int imag, const int id_thread){
    cdouble Hy_du;
    cdouble Hy_d;
    cdouble Hy_dl;
    cdouble dt = imag==0 ? cdouble(1.0,0.0)*_param->dt: cdouble(0.0,-1.0)*_param->dt_ITP;
    
    cdouble afield_j = (*Afield_j)[ti];
    cdouble bfield_k = (*Bfield_k)[ti];

    cdouble a = 1.0/(_dj*_dj)+0.5*afield_j*afield_j/(C*C);
    cdouble b = 1.0/(2.0*_dj*_dj);
    Hy_du = -b - I*1.0/(2.0*C*_dj)*afield_j-I/(2.0*_dj)*bfield_k*_i[i];
    Hy_dl = -b + I*1.0/(2.0*C*_dj)*afield_j+I/(2.0*_dj)*bfield_k*_i[i];
    for(int j=0;j<_nj;j++){
        Hy_d  =  a + 1.0/3.0*_potential_fn(_i[i],_j[j],_k[k],_t[ti]) + 1.0/3.0*0.125*bfield_k*bfield_k*(_i[i]*_i[i]+_j[j]*_j[j]);

        _Mj_du[id_thread*_nj + j] = I*Hy_du*dt/2.0;
        _Mj_d[id_thread*_nj + j]  = 1.0 + I*Hy_d*dt/2.0;
        _Mj_dl[id_thread*_nj + j] = I*Hy_dl*dt/2.0;
        _Mpj_du[id_thread*_nj + j] = -I*Hy_du*dt/2.0;
        _Mpj_d[id_thread*_nj + j]  = 1.0 - I*Hy_d*dt/2.0;
        _Mpj_dl[id_thread*_nj + j] = -I*Hy_dl*dt/2.0;
    }
    tridot(&_Mpj_dl[id_thread*_nj],
           &_Mpj_d[id_thread*_nj],
           &_Mpj_du[id_thread*_nj],
           psi_j_row, &_lhs_j[id_thread*_nj], _nj); 
    tdma(&_Mj_dl[id_thread*_nj],
         &_Mj_d[id_thread*_nj],
         &_Mj_du[id_thread*_nj],
         &_lhs_j[id_thread*_nj],&_res_j[id_thread*_nj],_nj);
    for(int j=0; j<_nj; j++){
        psi_j_row[j] = _res_j[id_thread*_nj + j];
    }
}
void Hamiltonian::step_k_XYZ(cdouble *psi_k_row, const int i, const int j, const int ti, const int imag, const int id_thread){
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
        Hz_d  =  a + 1.0/3.0*_potential_fn(_i[i],_j[j],_k[k],_t[ti]) + 1.0/3.0*0.125*bfield_k*bfield_k*(_i[i]*_i[i]+_j[j]*_j[j]);

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
           psi_k_row, &_lhs_k[id_thread*_nk], _nk); 
    tdma(&_Mk_dl[id_thread*_nk],
         &_Mk_d[id_thread*_nk],
         &_Mk_du[id_thread*_nk],
         &_lhs_k[id_thread*_nk],&_res_k[id_thread*_nk],_nk);
    for(int k=0; k<_nk; k++){
        psi_k_row[k] = _res_k[id_thread*_nk + k];
    }
}

cdouble Hamiltonian::ener_XYZ(cdouble ***psi){
    cdouble integral=0.0;

    cdouble *temp_x, *temp_y, *temp_z;
    cdouble *i_row, *j_row, *k_row;
    cdouble *Hx_du, *Hx_dl, *Hx_d;
    cdouble *Hz_du, *Hz_dl, *Hz_d; 
    cdouble *Hy_du, *Hy_dl, *Hy_d;

    double integral_z, integral_y, integral_x;
    // Apply Hz to wavefunction
    temp_z = new cdouble[_param->n_threads*_nk];
    k_row = new cdouble[_param->n_threads*_nk];
    Hz_dl = new cdouble[_param->n_threads*_nk];
    Hz_d = new cdouble[_param->n_threads*_nk];
    Hz_du = new cdouble[_param->n_threads*_nk];
    integral_z = 0.0;
    {    
    cdouble a = 1.0/(_dk*_dk);
    cdouble b = 1.0/(2.0*_dk*_dk);
    #pragma omp parallel for schedule(dynamic) collapse(1)  reduction(+: integral_z)
    for(int i=0;i<_ni;i++){
        for(int j=0;j<_nj;j++){
            int id = omp_get_thread_num();
            for(int k=0;k<_nk;k++){
                Hz_du[id*_nk + k] = -b;
                Hz_d[id*_nk + k]  =  a + 1.0/3.0*_potential_fn(_i[i],_j[j],_k[k],0);
                Hz_dl[id*_nk + k] = -b;
                k_row[id*_nk + k] = psi[i][j][k];
            }
            tridot(&Hz_dl[id*_nk], &Hz_d[id*_nk], &Hz_du[id*_nk], &k_row[id*_nk], &temp_z[id*_nk], _nk);
            for(int  k=0;k<_nk;k++){
                integral_z += (conj(psi[i][j][k])*temp_z[id*_nk + k]).real();
            }
        }
    }
    }
    delete[] temp_z;
    delete[] k_row;
    delete[] Hz_du;
    delete[] Hz_d;
    delete[] Hz_dl; 

    // Apply Hx to wavefunction
    temp_x = new cdouble[_param->n_threads*_ni];
    i_row = new cdouble[_param->n_threads*_ni];
    Hx_dl = new cdouble[_param->n_threads*_ni];
    Hx_d = new cdouble[_param->n_threads*_ni];
    Hx_du = new cdouble[_param->n_threads*_ni];
    integral_x = 0.0;
    {
    cdouble a = 1.0/(_di*_di);
    cdouble b = 1.0/(2.0*_di*_di); 
    #pragma omp parallel for schedule(dynamic) collapse(1)  reduction(+: integral_x)
    for(int j=0;j<_nj;j++){
        for(int k=0;k<_nk;k++){
            int id = omp_get_thread_num();
            for(int i=0;i<_ni;i++){
                Hx_du[id*_ni + i] = -b;
                Hx_d[id*_ni + i]  =  a + 1.0/3.0*_potential_fn(_i[i],_j[j],_k[k],0);
                Hx_dl[id*_ni + i] = -b;
                i_row[id*_ni + i] = psi[i][j][k];
            }
            tridot(&Hx_dl[id*_ni], &Hx_d[id*_ni], &Hx_du[id*_ni], &i_row[id*_ni], &temp_x[id*_ni], _ni);

            for(int i=0; i<_ni; i++){
                integral_x += (conj(psi[i][j][k])*temp_x[id*_ni + i]).real();
            }
        }
    }
    }
    delete[] temp_x;
    delete[] i_row;
    delete[] Hx_du;
    delete[] Hx_d;
    delete[] Hx_dl;

    // Apply Hy to wavefunction
    temp_y = new cdouble[_param->n_threads*_nj];
    j_row = new cdouble[_param->n_threads*_nj];
    Hy_dl = new cdouble[_param->n_threads*_nj];
    Hy_d = new cdouble[_param->n_threads*_nj];
    Hy_du = new cdouble[_param->n_threads*_nj];
    integral_y = 0.0;
    {
    cdouble a = 1.0/(_dj*_dj);
    cdouble b = 1.0/(2.0*_dj*_dj); 
    #pragma omp parallel for schedule(dynamic) collapse(1) reduction(+: integral_y)
    for(int i=0;i<_ni;i++){
        for(int k=0;k<_nk;k++){
            int id = omp_get_thread_num();
            for(int j=0;j<_nj;j++){
                Hy_du[id*_nj + j] = -b;
                Hy_d[id*_nj + j]  =  a + 1.0/3.0*_potential_fn(_i[i],_j[j],_k[k],0);
                Hy_dl[id*_nj + j] = -b;
                j_row[id*_nj + j] = psi[i][j][k];
            }
            tridot(&Hy_dl[id*_nj], &Hy_d[id*_nj], &Hy_du[id*_nj], &j_row[id*_nj], &temp_y[id*_nj], _nj);

            for(int j=0; j<_nj; j++){
                integral_y += (conj(psi[i][j][k])*temp_y[id*_nj + j]).real();
            }
        }
    }
    }
    delete[] temp_y;
    delete[] j_row;
    delete[] Hy_dl;
    delete[] Hy_d;
    delete[] Hy_du;
    // Integrate 
    integral = integral_x + integral_y + integral_z;
    integral *= _di*_dj*_dk;
    return integral; 
}

