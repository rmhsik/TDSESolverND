#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "debug.h"
#include "wavefunction.h"

WF::WF(){
        
}

WF::WF(Parameters param){
    _param = param;
}

void WF::set_geometry( double *i, double *k, const double di, const double dk){    
    _ni = _param.ni;
    _nk = _param.nk;

    _wf = new cdouble[_ni*_nk];
    _row = new cdouble[_ni];
    _col = new cdouble[_nk];

    for(int i=0; i<_ni; i++){
        for(int j=0; j<_nk;j++){
            _wf[i*_nk + j] = cdouble(0.0,0.0);
        }
    }

    _i = i; _k = k; _di = di; _dk = dk;
    
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
            _wf[i*_nk + j] = exp(-(_i[i]-i0)/sigma - (_k[j]-k0)/sigma);
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
    switch(_param.geometry){
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
    for(int i=0; i<_ni; i++){
        _wf[i*_nk + 0] = _wf[i*_nk + 0]*imask[i];
    }
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

cdouble WF::dipole(){
    cdouble sum = 0.0;
    switch(_param.geometry){
        case X:
            for(int i=0; i<_ni; i++){
                sum += conj(_wf[i*_nk + 0])*_i[i]*_wf[i*_nk + 0]*_di;
            }
            break;

        case XZ:
            for(int i=0; i<_ni; i++){
                for(int j=0; j<_nk; j++){
                    sum += conj(_wf[i*_nk + j])*_k[j]*_wf[i*_nk + j]*_di*_dk;
                }
            }
        case RZ:
            for(int i=0; i<_ni; i++){
                for(int j=0; j<_nk;j++){
                    sum += 2*M_PI*_i[i]*conj(_wf[i*_nk + j])*_k[j]*_wf[i*_nk+j]*_di*_dk;
                }
            }
    }   

    return sum; 
}

cdouble WF::acc(cdouble *dV){
    cdouble sum = cdouble(0.0,0.0);
    switch(_param.geometry){
        case X:
            for(int i=0; i<_ni; i++){
                sum += conj(_wf[i*_nk + 0])*(-1.0*dV[i*_nk + 0])*_wf[i*_nk + 0]*_di;
            }
            break;

        case XZ:
            for(int i=0; i<_ni; i++){
                for(int j=0; j<_nk; j++){
                    sum += conj(_wf[i*_nk + j])*(-1.0*dV[i*_nk + j])*_wf[i*_nk + j]*_di*_dk;
                }
            }
        case RZ:
            for(int i=0; i<_ni; i++){
                for(int j=0; j<_nk;j++){
                    sum += 2*M_PI*_i[i]*conj(_wf[i*_nk + j])*(-1.0*dV[i*_nk + j])*_wf[i*_nk + j]*_di*_dk;
                }
            }
    }   

    return sum; 
}
