#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include "debug.h"
#include "wavefunction.h"

WF::WF(){
        
}

WF::WF(Parameters *param){
    _param = param;
}

void WF::set_geometry( double *i, double *k, const double di, const double dk){    
    _ni = _param->ni;
    _nk = _param->nk;

    _wf = new cdouble[_ni*_nk];
    _row = new cdouble[_ni];
    _col = new cdouble[_nk];
    _wf_buf = new cdouble[_param->nt_diag*_ni*_nk];
    _diag_buf = new cdouble[_param->nt_diag];
    for(int i=0; i<_ni; i++){
        for(int j=0; j<_nk;j++){
            _wf[i*_nk + j] = cdouble(0.0,0.0);
        }
    }

    _i = i; _k = k; _di = di; _dk = dk;

    switch(_param->geometry){
        case X:
            _apply_mask = &WF::apply_mask_X;
            _pop = &WF::_pop_X;
            _pop_buf = &WF::_pop_buf_X;

            _acc_i = &WF::_acc_i_X;
            _acc_i_buf = &WF::_acc_i_buf_X;
            _acc_k = &WF::_acc_k_X;
            _acc_k_buf = &WF::_acc_k_buf_X;
            
            _dip_i = &WF::_dip_i_X;
            _dip_i_buf = &WF::_dip_i_buf_X;
            _dip_k = &WF::_dip_k_X;
            _dip_k_buf = &WF::_dip_k_buf_X;
            break;

        case XZ:
            _apply_mask = &WF::apply_mask_XZ;
            _pop = &WF::_pop_XZ;
            _pop_buf = &WF::_pop_buf_XZ;
        
            _acc_i = &WF::_acc_i_XZ;
            _acc_i_buf = &WF::_acc_i_buf_XZ;
            _acc_k = &WF::_acc_k_XZ;
            _acc_k_buf = &WF::_acc_k_buf_XZ;

            _dip_i = &WF::_dip_i_XZ;
            _dip_i_buf = &WF::_dip_i_buf_XZ;
            _dip_k = &WF::_dip_k_XZ;
            _dip_k_buf = &WF::_dip_k_buf_XZ;
            break;
        case RZ:
            _apply_mask = &WF::apply_mask_RZ;
            _pop = &WF::_pop_RZ;
            _pop_buf = &WF::_pop_buf_RZ;

            _acc_i = &WF::_acc_i_RZ;
            _acc_i_buf = &WF::_acc_i_buf_RZ;
            _acc_k = &WF::_acc_k_RZ;
            _acc_k_buf = &WF::_acc_k_buf_RZ;

            _dip_i = &WF::_dip_i_RZ;
            _dip_i_buf = &WF::_dip_i_buf_RZ;
            _dip_k = &WF::_dip_k_RZ;
            _dip_k_buf = &WF::_dip_k_buf_RZ;
            break;
    }
}

void WF::gaussian(double i0, double k0, double sigma){
    for(int i=0; i<_ni; i++){
        for(int j=0; j<_nk;j++){
            _wf[i*_nk + j] = exp(-(_i[i]-i0)*(_i[i]-i0)/sigma - (_k[j]-k0)*(_k[j]-k0)/sigma);
        }
    }

}

void WF::exponential(double i0, double k0, double sigma){
    for(int i=0; i<_ni; i++){
        for(int j=0; j<_nk;j++){
            _wf[i*_nk + j] = exp(-sqrt((_i[i]-i0)*(_i[i]-i0) + (_k[j]-k0)*(_k[j]-k0)));
        }
    }
}


void WF::save_wf(std::string path){
    std::ofstream outfile;
    outfile.open(path);
    if(outfile.is_open()){
        for(int i=0;i<_ni;i++){
            for(int j=0; j<_nk; j++){
                std::stringstream str;
                str<<std::fixed<<std::setprecision(12);
                str<<_wf[i*_nk + j];
                outfile<<str.str();
            }
            outfile<<std::endl;
        }
    }
    else{
        debug1("[WF->save_wf] Error opening output file");
    }    
}

void WF::save_wf2(std::string path){
    std::ofstream outfile;
    outfile.open(path);
    if(outfile.is_open()){
        for(int i=0;i<_ni;i++){
            for(int j=0; j<_nk; j++){
                std::stringstream str;
                str<<std::fixed<<std::setprecision(12);
                str<<(double)abs(_wf[i*_nk + j]*_wf[i*_nk + j]);
                outfile<<str.str();
            }
            outfile<<std::endl;
        }
    }
    else{
        debug1("[WF->save_wf2] Error opening output file");
    }    
}

cdouble WF::norm(){
    cdouble integral = 0.0;
    switch(_param->geometry){
        case X:
            for(int i=0; i<_ni;i++){
                integral += _wf[i*_nk + 0]*conj(_wf[i*_nk + 0])*_di;
            }
            break;
        case XZ:
            for(int i=0; i<_ni;i++){
                for(int j=0;j<_nk;j++){
                    integral += _wf[i*_nk + j]*conj(_wf[i*_nk + j])*_di*_dk;
                }
            }
            break;
        case RZ:
            for(int i=0; i<_ni;i++){
                for(int j=0;j<_nk;j++){
                    integral += 2*M_PI*_i[i]*_wf[i*_nk + j]*conj(_wf[i*_nk + j])*_di*_dk;
                }
            }
            break;
    }
    return sqrt(integral);
}

void WF::apply_mask(cdouble *imask, cdouble *kmask){
    (this->*(this->_apply_mask))(imask, kmask);
}


cdouble* WF::get(){
    return _wf;
}

cdouble* WF::row(int k){
    for(int i=0; i<_ni; i++){
        _row[i] = _wf[i*_nk + k];
    }
    return _row;
}

cdouble* WF::col(int i){
    for(int k=0; k<_nk; k++){
        _col[k] = _wf[i*_nk + k];
    }
    return _col;
}

void WF::set_row(cdouble* row, int k){
    for(int i=0; i<_ni; i++){
        _wf[i*_nk + k] = row[i];
    }
}

void WF::set_col(cdouble* col, int i){
    for(int k=0; k<_nk;k++){
        _wf[i*_nk + k] = col[k];
    }
}

void WF::set(cdouble* arr){
    for(int i=0; i<_ni;i++){
        for(int j=0;j<_nk;j++){
            _wf[i*_nk + j] = arr[i*_nk + j];
        }
    }
}

void WF::set_to_buf(const int idx){
    std::memcpy(&_wf_buf[idx*_ni*_nk],_wf,_ni*_nk*sizeof(cdouble));
}

void WF::set_col_buf(cdouble* col, const int i, const int idx){
    for(int k=0; k<_nk; k++)
        _wf_buf[idx*_ni*_nk + i*_nk + k] = col[k];
}

void WF::set_row_buf(cdouble *row, const int k, const int idx){
    for(int i=0;i<_ni;i++)
        _wf_buf[idx*_ni*_nk + i*_nk +k] = row[i];
}

void WF::set_col_buf_mask(cdouble* col, cdouble* kmask, const int i, const int idx){
    for(int k=0; k<_nk; k++)
        _wf_buf[idx*_ni*_nk + i*_nk + k] = col[k]*kmask[k];
}

void WF::set_row_buf_mask(cdouble* row, cdouble* imask, const int k, const int idx){
    for(int i=0; i<_ni; i++)
        _wf_buf[idx*_ni*_nk + i*_nk + k] = row[i]*imask[i];
}


void WF::get_from_buf(cdouble* arr, const int idx){
    std::memcpy(arr, &_wf_buf[idx*_ni*_nk],_ni*_nk*sizeof(cdouble));
}

cdouble* WF::get_buf(){
    return _wf_buf;
}

cdouble* WF::get_diag_buf(){
    return _diag_buf;
}

void WF::set_dpotential_i(cdouble* dV){
    _dV_i = dV;
}

void WF::set_dpotential_k(cdouble* dV){
    _dV_k = dV;
}

cdouble WF::operator()(int i, int j){
    return _wf[i*_nk + j];
}

void WF::operator/=(cdouble val){
    for(int i=0; i<_ni;i++){
        for(int j=0; j<_nk;j++){
            _wf[i*_nk + j] /= val;
        }
    }
}

cdouble WF::_dip_0(){return 0.0;}
void WF::_dip_buf_0(){}
cdouble WF::_acc_0(){return 0.0;}
void WF::_acc_buf_0(){}
cdouble WF::_pop_0(double imin, double imax, double kmin, double kmax){return 0.0;}
void WF::_pop_buf_0(double imin, double imax, double kmin, double kmax){}

cdouble WF::acc_i(){
    cdouble sum = (this->*(this->_acc_i))();
    return sum; 
}

cdouble WF::acc_k(){
    cdouble sum = (this->*(this->_acc_k))();
    return sum;
}

cdouble WF::dip_i(){
    cdouble sum = (this->*(this->_dip_i))();
    return sum;
}  

cdouble WF::dip_k(){
    cdouble sum = (this->*(this->_dip_k))();
    return sum;
}

void WF::acc_i_buf(){
    (this->*(this->_acc_i_buf))();
}
void WF::acc_k_buf(){
    (this->*(this->_acc_k_buf))();
}

void WF::dip_i_buf(){
    (this->*(this->_dip_i_buf))();
}
void WF::dip_k_buf(){
    (this->*(this->_dip_k_buf))();
}

cdouble WF::pop(double imin, double imax, double kmin, double kmax){
    return (this->*(this->_pop))(imin, imax, kmin, kmax);
}

void WF::pop_buf(double imin, double imax, double kmin, double kmax){
    (this->*(this->_pop_buf))(imin, imax, kmin, kmax);
}
