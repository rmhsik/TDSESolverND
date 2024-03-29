#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <vector>
#include "debug.h"
#include "utils.h"
#include "wavefunction.h"

WF::WF(){
        
}

WF::WF(Parameters *param){
    _param = param;
}

void WF::set_geometry( double *i,double *k, const double di, const double dk){    
    _ni = _param->ni;
    _nk = _param->nk;
    
    _wf = alloc2d<cdouble>(_ni, _nk);
    _wf_buf = alloc3d<cdouble>(_param->nt_diag,_ni, _nk);

    _i_row = new cdouble[_ni];
    _k_row = new cdouble[_nk];

    _diag_buf = new cdouble[_param->nt_diag];

    for(int i=0; i<_ni; i++){
        for(int k=0;k<_nk;k++){
            _wf[i][k] = 0.0;
        }
    }

    _i = i; _k = k; _di = di; _dk = dk;

    switch(_param->geometry){
        case XZ:
            _geom_XZ();
            break;
    }
}

void WF::gaussian(double i0, double k0, double sigma){
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<_ni; i++){
    	for(int k=0; k<_nk; k++){
        	_wf[i][k] = exp(-(_i[i]-i0)*(_i[i]-i0)/sigma - (_k[k]-k0)*(_k[k]-k0)/sigma);
        }
    }
}

void WF::gaussian_anti(double i0, double k0, double sigma){
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<_ni; i++){
		for(int k=0; k<_nk; k++){
     		_wf[i][k] = _k[k]*exp(-(_i[i]-i0)*(_i[i]-i0)/sigma- (_k[k]-k0)*(_k[k]-k0)/sigma);
        }
    }
}

void WF::exponential(double i0, double k0, double sigma){
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<_ni; i++){
        for(int k=0; k<_nk;k++){
            _wf[i][k] = exp(-sqrt((_i[i]-i0)*(_i[i]-i0) + (_k[k]-k0)*(_k[k]-k0)));
        }
    }
}


cdouble WF::norm(){
    cdouble integral = 0.0;
    switch(_param->geometry){
        case XZ:
            for(int i=0; i<_ni;i++){
                for(int k=0;k<_nk;k++){
                    integral += _wf[i][k]*conj(_wf[i][k])*_di*_dk;
                }
            }
            break;
    }
    return sqrt(integral);
}

cdouble WF::norm_buf(const int idx){
    cdouble integral = 0.0;
    switch(_param->geometry){
        case XZ:
            for(int i=0; i<_ni;i++){
                for(int k=0;k<_nk;k++){
                    integral += _wf_buf[idx][i][k]*conj(_wf_buf[idx][i][k])*_di*_dk;
                }
            }
            break;
    }
    return sqrt(integral);
}

void WF::apply_mask(cdouble *imask, cdouble *kmask){
    (this->*(this->_apply_mask))(imask, kmask);
}


cdouble** WF::get(){
    return _wf;
}

cdouble* WF::i_row(int k){
    for(int i=0; i<_ni; i++){
        _i_row[i] = _wf[i][k];
    }
    return _i_row;
}

cdouble* WF::k_row(int i){
    for(int k=0; k<_nk;k++){
        _k_row[k] = _wf[i][k];
    }
    return _k_row;
}

void WF::set_i_row(cdouble* i_row, int k){
    for(int i=0; i<_ni; i++){
        _wf[i][k] = i_row[i];
    }
}

void WF::set_i_row_mask(cdouble* i_row, cdouble* imask, int k){
    for(int i=0; i<_ni; i++){
        _wf[i][k] = i_row[i]*imask[i];
    }
}

void WF::set_k_row(cdouble* k_row, int i){
    for(int k=0; k<_nk; k++){
        _wf[i][k] = k_row[k];
    }
}

void WF::set_k_row_mask(cdouble* k_row, cdouble *kmask, int i){
    for(int k=0; k<_nk; k++){
        _wf[i][k] = k_row[k]*kmask[k];
    }
}

void WF::get_i_row(cdouble* i_row,int k){
    for(int i=0; i<_ni; i++){
        i_row[i] = _wf[i][k];
    }
}

void WF::get_k_row(cdouble* k_row, int i){
    for(int k=0; k<_nk; k++){
        k_row[k] = _wf[i][k];
    }
}
void WF::set(cdouble** arr){
    for(int i=0; i<_ni;i++){
        for(int k=0; k<_nk;k++){
            _wf[i][k] = arr[i][k];
        }
    }
}

void WF::set_to_buf(const int idx){
    for(int i=0;i<_ni;i++){
        for(int k=0;k<_nk;k++){
            _wf_buf[idx][i][k] = _wf[i][k];
        }
    }
    //std::memcpy(_wf_buf[idx],_wf,_ni*_nk*sizeof(cdouble));
}

void WF::set_i_row_buf(cdouble* i_row, const int k, const int idx){
    for(int i=0; i<_ni; i++){
        _wf_buf[idx][i][k] = i_row[i];
    }
}

void WF::set_k_row_buf(cdouble *k_row, const int i, const int idx){
    for(int k=0; k<_nk;k++)
        _wf_buf[idx][i][k] = k_row[k];
}

void WF::get_i_row_buf(cdouble* i_row, const int k, const int idx){
    for(int i=0; i<_ni; i++){
        i_row[i] = _wf_buf[idx][i][k];
    }
}

void WF::get_k_row_buf(cdouble *k_row, const int i, const int idx){
    for(int k=0; k<_nk;k++)
        k_row[k] = _wf_buf[idx][i][k];
}


void WF::set_i_row_buf_mask(cdouble* i_row, cdouble* imask,const int k, const int idx){
    for(int i=0; i<_ni; i++)
        _wf_buf[idx][i][k] = i_row[i]*imask[i];
}

void WF::set_k_row_buf_mask(cdouble* k_row, cdouble* kmask, const int i, const int idx){
    for(int k=0; k<_nk; k++)
        _wf_buf[idx][i][k] = k_row[k]*kmask[k];
}



void WF::get_from_buf(cdouble** arr, const int idx){
    for(int i=0; i<_ni; i++){
        for(int k=0;k<_nk; k++){
            arr[i][k] = _wf_buf[idx][i][k];
        }
    }
}

void WF::anti_sym_k(){
    for(int i=0; i<_ni; i++){
        for(int k=0;k<_nk/2;k++){
            _wf[i][k] = (_wf[i][k] - _wf[i][_nk-k-1])/2.0;           
            _wf[i][_nk-k-1] = -_wf[i][k];
        }
    }
}

cdouble*** WF::get_buf(){
    return _wf_buf;
}

cdouble WF::operator()(int i, int k){
    return _wf[i][k];
}

void WF::operator/=(cdouble val){
    #pragma omp parallel for collapse(2) schedule(static)
    for(int i=0; i<_ni;i++){
        for(int k=0; k<_nk;k++){
            _wf[i][k] /= val;
        }
    }
}

void WF::save_wf2(std::string name){
}

WF::~WF(){
    free2d(&_wf,_ni,_nk);
    free3d(&_wf_buf, _param->nt_diag, _ni, _nk);
    //if(_param->geometry != XYZ)
    //    free3d(&_wf_0,_ni,_nk,_nk);
    //else
    //    free3d(&_wf_0,1,1,1);
    delete[] _diag_buf;
    delete[] _i_row;
    delete[] _k_row;
}
