#include <cmath>
#include "diagnostics.h"
#include "utils.h"

Diagnostics::Diagnostics(){

}

Diagnostics::Diagnostics(Parameters *param){
   _param = param; 
   _acc_i = new cdouble[_param->nt];
   _acc_k = new cdouble[_param->nt];
   _dip_i = new cdouble[_param->nt];
   _dip_k = new cdouble[_param->nt];
   _pop   = new cdouble[_param->nt];
}

void Diagnostics::set_geometry(double *i, double *k, const double di, const double dk){
    _i = i; _k = k; _di = di; _dk = dk;
   _ni = _param->ni;
   _nk = _param->nk;
}

cdouble* Diagnostics::get_acc_i(){
    return _acc_i;
}

cdouble* Diagnostics::get_acc_k(){
    return _acc_k;
}

cdouble* Diagnostics::get_dip_i(){
    return _dip_i;
}

cdouble* Diagnostics::get_dip_k(){
    return _dip_k;
}

cdouble* Diagnostics::get_pop(){
    return _pop;
}

void Diagnostics::write_diagnostics(){
    write_array(_acc_i, _param->nt, _param->acc_i_path);
    write_array(_acc_k, _param->nt, _param->acc_k_path);
    write_array(_dip_i, _param->nt, _param->dip_i_path);
    write_array(_dip_k, _param->nt, _param->dip_k_path);
    write_array(_pop,   _param->nt, _param->pop_path);
}


