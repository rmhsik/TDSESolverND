#include <iostream>
#include "hamiltonian.h"

void Hamiltonian::_allocate_X(){
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

}

void Hamiltonian::step_i_X(cdouble *psi, const int j, const int ti, const int imag,const int id_thread){
    cdouble H_du;
    cdouble H_d;
    cdouble H_dl;
    cdouble dt = imag==0 ? cdouble(1.0,0.0)*_param->dt: cdouble(0.0,-1.0)*_param->dt_ITP;

    cdouble afield_i = (*Afield_i)[ti];

    cdouble a = 1.0/(_di*_di);
    cdouble b = 1.0/(2.0*_di*_di); 
    for(int i=0;i<_ni;i++){
        H_du = -b + I*1.0/(2.0*C*_di)*afield_i;
        H_d  =  a + _potential_fn(_i[i],0,_t[ti]) + 0.5*afield_i*afield_i/(C*C);
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
        H_d[i]  =  a + _potential_fn(_i[i],0,0);
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


