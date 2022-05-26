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

void WF::set_geometry( double *i, double *j, double *k, const double di, const double dj, const double dk){    
    _ni = _param->ni;
    _nj = _param->nj;
    _nk = _param->nk;
    
    _wf = alloc3d<cdouble>(_ni, _nj, _nk);
    if(_param->geometry == XYZ)
        _wf_buf = alloc4d<cdouble>(_ni, _nj, _nk, 1);
    else
        _wf_buf = alloc4d<cdouble>(_ni, _nj, _nk, _param->nt_diag);

    _eigen_wf = alloc4d<cdouble>(_ni,_nj,_nk,5);

    _i_row = new cdouble[_ni];
    _j_row = new cdouble[_nj];
    _k_row = new cdouble[_nk];

    _diag_buf = new cdouble[_param->nt_diag];

    for(int i=0; i<_ni; i++){
        for(int j=0; j<_nj;j++){
            for(int k=0;k<_nk;k++){
                _wf[i][j][k] = 0.0;
            }
        }
    }


    for(int n=0;n<5;n++){
        for(int i=0;i<_ni;i++){
            for(int j=0;j<_nj;j++){
                for(int k=0;k<_nk;k++){
                    //std::cout<<i<<" "<<k<<std::endl;
                    _eigen_wf[n][i][j][k] = 0.0;
                }
            }
        }
    }

    _i = i; _j = j; _k = k; _di = di; _dj = dj; _dk = dk;

    switch(_param->geometry){
        case X:
            _geom_X();
            break;

        case XZ:
            _geom_XZ();
            break;
        case RZ:
            _geom_RZ();
            break;
        case XYZ:
            _geom_XYZ();
            break;
    }
}

void WF::gaussian(double i0, double j0, double k0, double sigma){
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<_ni; i++){
        for(int j=0; j<_nj;j++){
            for(int k=0; k<_nk; k++){
                _wf[i][j][k] = exp(-(_i[i]-i0)*(_i[i]-i0)/sigma- (_j[j]-j0)*(_j[j]-j0)/sigma - (_k[k]-k0)*(_k[k]-k0)/sigma);
            }
        }
    }
}

void WF::gaussian_anti(double i0, double j0, double k0, double sigma){
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<_ni; i++){
        for(int j=0; j<_nj;j++){
            for(int k=0; k<_nk; k++){
                _wf[i][j][k] = _k[k]*exp(-(_i[i]-i0)*(_i[i]-i0)/sigma- (_j[j]-j0)*(_j[j]-j0)/sigma - (_k[k]-k0)*(_k[k]-k0)/sigma);
            }
        }
    }
}

void WF::exponential(double i0, double j0, double k0, double sigma){
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<_ni; i++){
        for(int j=0; j<_nj;j++){
            for(int k=0; k<_nk;k++){
                _wf[i][j][k] = exp(-sqrt((_i[i]-i0)*(_i[i]-i0) + (_j[j]-j0)*(_j[j]-j0) + (_k[k]-k0)*(_k[k]-k0)));
            }
        }
    }
}

cdouble WF::norm(){
    cdouble integral = 0.0;
    switch(_param->geometry){
        case X:
            for(int i=0; i<_ni;i++){
                integral += _wf[i][0][0]*conj(_wf[i][0][0])*_di;
            }
            break;
        case XZ:
            for(int i=0; i<_ni;i++){
                for(int k=0;k<_nk;k++){
                    integral += _wf[i][0][k]*conj(_wf[i][0][k])*_di*_dk;
                }
            }
            break;
        case RZ:
            for(int i=0; i<_ni;i++){
                for(int k=0;k<_nk;k++){
                    integral += 2*M_PI*_i[i]*_wf[i][0][k]*conj(_wf[i][0][k])*_di*_dk;
                }
            }
            break;
        case XYZ:
            double sum = 0.0;
            //#pragma omp parallel for reduction(+:sum)
            for(int i=0;i<_ni;i++){
                for(int j=0;j<_nj;j++){
                    for(int k=0;k<_nk;k++){
                        sum += (_wf[i][j][k]*conj(_wf[i][j][k])).real()*_di*_dk*_dj;
                    }
                }
            }
            integral = sum;
            break;
    }
    return sqrt(integral);
}

void WF::apply_mask(cdouble *imask, cdouble *jmask, cdouble *kmask){
    (this->*(this->_apply_mask))(imask, jmask, kmask);
}


cdouble*** WF::get(){
    return _wf;
}

cdouble* WF::i_row(int j, int k){
    for(int i=0; i<_ni; i++){
        _i_row[i] = _wf[i][j][k];
    }
    return _i_row;
}

cdouble* WF::j_row(int i, int k){
    for(int j=0; j<_nj; j++){
        _j_row[j] = _wf[i][j][k];
    }
    return _j_row;
}

cdouble* WF::k_row(int i, int j){
    for(int k=0; k<_nk;k++){
        _k_row[k] = _wf[i][j][k];
    }
    return _k_row;
}

void WF::set_i_row(cdouble* i_row, int j, int k){
    for(int i=0; i<_ni; i++){
        _wf[i][j][k] = i_row[i];
    }
}

void WF::set_i_row_mask(cdouble* i_row, cdouble* imask, int j, int k){
    for(int i=0; i<_ni; i++){
        _wf[i][j][k] = i_row[i]*imask[i];
    }
}

void WF::set_j_row(cdouble* j_row, int i, int k){
    for(int j=0; j<_nj;j++){
        _wf[i][j][k] = j_row[j];
    }
}

void WF::set_j_row_mask(cdouble* j_row, cdouble *jmask, int i, int k){
    for(int j=0; j<_nj;j++){
        _wf[i][j][k] = j_row[j]*jmask[j];
    }
}

void WF::set_k_row(cdouble* k_row, int i, int j){
    for(int k=0; k<_nk; k++){
        _wf[i][j][k] = k_row[k];
    }
}

void WF::set_k_row_mask(cdouble* k_row, cdouble *kmask, int i, int j){
    for(int k=0; k<_nk; k++){
        _wf[i][j][k] = k_row[k]*kmask[k];
    }
}

void WF::get_i_row(cdouble* i_row, int j, int k){
    for(int i=0; i<_ni; i++){
        i_row[i] = _wf[i][j][k];
    }
}

void WF::get_j_row(cdouble* j_row, int i, int k){
    for(int j=0; j<_nj;j++){
        j_row[j] = _wf[i][j][k];
    }
}

void WF::get_k_row(cdouble* k_row, int i, int j){
    for(int k=0; k<_nk; k++){
        k_row[k] = _wf[i][j][k];
    }
}
void WF::set(cdouble*** arr){
    for(int i=0; i<_ni;i++){
        for(int j=0;j<_nj;j++){
            for(int k=0; k<_nk;k++){
                _wf[i][j][k] = arr[i][j][k];
            }
        }
    }
}

void WF::set_to_buf(const int idx){
    for(int i=0;i<_ni;i++){
        for(int j=0;j<_nj;j++){
            for(int k=0;k<_nk;k++){
                _wf_buf[idx][i][j][k] = _wf[i][j][k];
            }
        }
    }
    //std::memcpy(_wf_buf[idx],_wf,_ni*_nk*sizeof(cdouble));
}

void WF::set_i_row_buf(cdouble* i_row, const int j, const int k, const int idx){
    for(int i=0; i<_ni; i++){
        _wf_buf[idx][i][j][k] = i_row[i];
    }
}

void WF::set_j_row_buf(cdouble *j_row, const int i, const int k, const int idx){
    for(int j=0;j<_nj;j++)
        _wf_buf[idx][i][j][k] = j_row[j];
}

void WF::set_k_row_buf(cdouble *k_row, const int i, const int j, const int idx){
    for(int k=0; k<_nk;k++)
        _wf_buf[idx][i][j][k] = k_row[k];
}

void WF::get_i_row_buf(cdouble* i_row, const int j, const int k, const int idx){
    for(int i=0; i<_ni; i++){
        i_row[i] = _wf_buf[idx][i][j][k];
    }
}

void WF::get_j_row_buf(cdouble *j_row, const int i, const int k, const int idx){
    for(int j=0;j<_nj;j++)
        j_row[j] = _wf_buf[idx][i][j][k];
}

void WF::get_k_row_buf(cdouble *k_row, const int i, const int j, const int idx){
    for(int k=0; k<_nk;k++)
        k_row[k] = _wf_buf[idx][i][j][k];
}


void WF::set_i_row_buf_mask(cdouble* i_row, cdouble* imask, const int j, const int k, const int idx){
    for(int i=0; i<_ni; i++)
        _wf_buf[idx][i][j][k] = i_row[i]*imask[i];
}

void WF::set_j_row_buf_mask(cdouble* j_row, cdouble* jmask, const int i, const int k, const int idx){
    for(int j=0; j<_nj; j++)
        _wf_buf[idx][i][j][k] = j_row[j]*jmask[j];
}

void WF::set_k_row_buf_mask(cdouble* k_row, cdouble* kmask, const int i, const int j, const int idx){
    for(int k=0; k<_nk; k++)
        _wf_buf[idx][i][j][k] = k_row[k]*kmask[k];
}



void WF::get_from_buf(cdouble*** arr, const int idx){
    for(int i=0; i<_ni; i++){
        for(int j=0; j<_nj; j++){
            for(int k=0;k<_nk; k++){
                arr[i][j][k] = _wf_buf[idx][i][j][k];
            }
        }
    }
}

void WF::anti_sym_k(){
    for(int i=0; i<_ni; i++){
        for(int j=0; j<_nj; j++){
            for(int k=0;k<_nk/2;k++){
                _wf[i][j][k] = (_wf[i][j][k] - _wf[i][j][_nk-k-1])/2.0;           
                _wf[i][j][_nk-k-1] = -_wf[i][j][k];
            }
        }
    }
}

cdouble**** WF::get_buf(){
    return _wf_buf;
}

cdouble* WF::get_diag_buf(){
    return _diag_buf;
}

void WF::set_to_eigen(int n){
    for(int i=0; i<_ni; i++ ){
        for(int j=0; j<_nj; j++){
            for(int k=0; k<_nk; k++){
                _eigen_wf[n][i][j][k] = _wf[i][j][k];
            }
        }
    }
}

cdouble*** WF::get_eigen(int n){
    return _eigen_wf[n];
}


void WF::set_to_ground(){
    for(int i=0;i<_ni;i++){
        for(int j=0;j<_nj;j++){
            for(int k=0;k<_nk;k++){
                _eigen_wf[0][i][j][k] = _wf[i][j][k];
            }
        }
    }
}

cdouble*** WF::get_ground(){
    return _eigen_wf[0];
}

cdouble WF::operator()(int i, int j, int k){
    return _wf[i][j][k];
}

cdouble WF::project(cdouble ***phi){
    cdouble integral = 0.0;
    switch(_param->geometry){
        case X:
            for(int i=0; i<_ni;i++){
                integral += _wf[i][0][0]*conj(phi[i][0][0])*_di;
            }
            break;
        case XZ:
            for(int i=0; i<_ni;i++){
                for(int k=0;k<_nk;k++){
                    integral += _wf[i][0][k]*conj(phi[i][0][k])*_di*_dk;
                }
            }
            break;
        case RZ:
            for(int i=0; i<_ni;i++){
                for(int k=0;k<_nk;k++){
                    integral += 2*M_PI*_i[i]*_wf[i][0][k]*conj(phi[i][0][k])*_di*_dk;
                }
            }
            break;
        case XYZ:
            //#pragma omp parallel for reduction(+:sum)
            for(int i=0;i<_ni;i++){
                for(int j=0;j<_nj;j++){
                    for(int k=0;k<_nk;k++){
                        integral += (_wf[i][j][k]*conj(phi[i][j][k]))*_di*_dk*_dj;
                    }
                }
            }
            break;
    }
    return integral;
}

void WF::grand_schmidt(){
    std::vector<cdouble> proj;
    for(int i=0;i<5;i++){
        proj.push_back(project(_eigen_wf[i]));
    }


    std::cout<<"Projection: "<<proj[0]<<" "<<proj[1]<<" "<<proj[2]<<" "<<proj[3]<<std::endl;
    for (int i=0; i<_ni; i++){
        for (int j=0; j<_nj; j++){
            for (int k=0; k<_nk; k++){
                for(int n=0; n<5; n++){
                    _wf[i][j][k] -= proj[n]*_eigen_wf[n][i][j][k];
                } 
            }
        }
    }
}


void WF::operator/=(cdouble val){
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for(int i=0; i<_ni;i++){
        for(int j=0; j<_nj;j++){
            for(int k=0; k<_nk;k++){
                _wf[i][j][k] /= val;
            }
        }
    }
}

void WF::save_wf2(std::string name){

    switch(_param->geometry){
        case XYZ:
            cdouble value;
            std::ofstream outfile_X;
            std::ofstream outfile_Y;
            std::ofstream outfile_Z;
            outfile_X.open("results/" + name + "_X.dat");
            outfile_Y.open("results/" + name + "_Y.dat");
            outfile_Z.open("results/" + name + "_Z.dat");

            for(int j=0;j<_nj;j++){
                for(int k=0; k<_nk;k++){
                    value = _wf[_ni/2][j][k]*conj(_wf[_ni/2][j][k]);
                    std::ostringstream doubleStr;
                    doubleStr<<std::fixed<<std::setprecision(12);
                    doubleStr<<std::real(value);
                    outfile_X<<doubleStr.str()<<std::endl;
                }
            }

            for(int i=0;i<_ni;i++){
                for(int k=0; k<_nk;k++){
                    value = _wf[i][_nj/2][k]*conj(_wf[i][_nj/2][k]);
                    std::ostringstream doubleStr;
                    doubleStr<<std::fixed<<std::setprecision(12);
                    doubleStr<<std::real(value);
                    outfile_Y<<doubleStr.str()<<std::endl;
                }
            }

            for(int i=0;i<_ni;i++){
                for(int j=0; j<_nj;j++){
                    value = _wf[i][j][_nk/2]*conj(_wf[i][j][_nk/2]);
                    std::ostringstream doubleStr;
                    doubleStr<<std::fixed<<std::setprecision(12);
                    doubleStr<<std::real(value);
                    outfile_Z<<doubleStr.str()<<std::endl;
                }
            }
            outfile_X.close();
            outfile_Y.close();
            outfile_Z.close();                 
            break;
    }
}

WF::~WF(){
    free3d(&_wf,_ni,_nj,_nk);
    if(_param->geometry == XYZ)
        free4d(&_wf_buf,_ni,_nj,_nk, 1);
    else
        free4d(&_wf_buf, _ni, _nj, _nk, _param->nt_diag);
    free4d(&_eigen_wf, 3, _ni, _nj, _nk);
    //if(_param->geometry != XYZ)
    //    free3d(&_wf_0,_ni,_nk,_nk);
    //else
    //    free3d(&_wf_0,1,1,1);
    delete[] _diag_buf;
    delete[] _i_row;
    delete[] _j_row;
    delete[] _k_row;
}
