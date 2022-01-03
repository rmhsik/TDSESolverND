#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "debug.h"
#include "wavefunction.h"

WF::WF(){
        
}

WF::WF(Parameters param,double *i, double *k, const double di, const double dk ){
    _param = param;
    
    _ni = _param.ni;
    if(_param.geometry == X){
        _nk=1;
    }
    else{
        _nk = _param.nk;
    }

    _wf = new cdouble*[_ni];
    for (int i=0;i<_ni;i++){
        _wf[i] = new cdouble[_nk];
    }

    for(int i=0; i<_ni; i++){
        for(int j=0; j<_nk;j++){
            _wf[i][j] = cdouble(0.0,0.0);
        }
    }

    _i = i; _k = k; _di = di; _dk = dk;

}

void WF::gaussian(double i0, double k0, double sigma){
    for(int i=0; i<_ni; i++){
        for(int j=0; j<_nk;j++){
            _wf[i][j] = exp(-(_i[i]-i0)*(_i[i]-i0)/sigma - (_k[j]-k0)*(_k[j]-k0)/sigma);
        }
    }

}

void WF::exponential(double i0, double k0, double sigma){
    //TODO
}


void WF::save_wf(std::string path){
    std::ofstream outfile;
    outfile.open(path);
    if(outfile.is_open()){
        for(int i=0;i<_ni;i++){
            for(int j=0; j<_nk; j++){
                std::stringstream str;
                str<<std::fixed<<std::setprecision(12);
                str<<_wf[i][j];
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
                str<<(double)abs(_wf[i][j]*_wf[i][j]);
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
                integral += _wf[i][0]*conj(_wf[i][0])*_di;
            }
            break;
        case XZ:
            for(int i=0; i<_ni;i++){
                for(int j=0;j<_nk;j++){
                    integral += _wf[i][j]*conj(_wf[i][j])*_di*_dk;
                }
            }
            break;
        case RZ:
            for(int i=0; i<_ni;i++){
                for(int j=0;j<_nk;j++){
                    integral += 2*M_PI*_i[i]*_wf[i][j]*conj(_wf[i][j])*_di*_dk;
                }
            }
            break;
    }
    return sqrt(integral);
}

cdouble** WF::get(){
    return _wf;
}

cdouble WF::operator()(int i, int j){
    return _wf[i][j];
}

void WF::operator/=(cdouble val){
    for(int i=0; i<_ni;i++){
        for(int j=0; j<_nk;j++){
            _wf[i][j] /= val;
        }
    }
}
