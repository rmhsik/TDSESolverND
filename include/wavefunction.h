#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H


#define cdouble std::complex<double>
#define I std::complex<double>(0.0,1.0)
#include <complex>
#include <string>
#include "parameters.h"

class WF{
    private:
        cdouble **_wf;
        Parameters _param;
        int _ni, _nk;
        double *_i, *_k, _di, _dk;
    public:
        WF();
        WF(Parameters param, double *i, double *k, const double di, const double _dk);
        void gaussian(double i0, double k0, double sigma);
        void exponential(double i0, double k0, double sigma);
        cdouble** get();
        cdouble norm();
        void save_wf(std::string path);
        void save_wf2(std::string path);
        cdouble operator()(int i, int j);
        void operator/= (cdouble val);        
};



#endif
