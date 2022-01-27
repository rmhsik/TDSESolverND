#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include "wavefunction.h"
#include "hamiltonian.h"
#include "parameters.h"

class Diagnostics{
    protected:
        Parameters *_param;
        cdouble *_acc_i, *_acc_k;
        cdouble *_dip_i, *_dip_k;
        cdouble *_pop;
     
        int _ni, _nk;
        double _di, _dk; 
        double *_i, *_k;

        Hamiltonian *_ham;
        WF *_wf;
        virtual void _acc_i_buf(const int idx){};
        virtual void _acc_k_buf(const int idx){};
        virtual void _dip_i_buf(const int idx){};
        virtual void _dip_k_buf(const int idx){};
        virtual void _pop_buf(const int idx){};
    public:
        Diagnostics();
        Diagnostics(Parameters * param);
        void set_geometry(double *i, double *k, const double di, const double dk);
        void set_ham(Hamiltonian *ham);
        void set_wf(WF *wf); 
        void run_diagnostics(const int idx);
        cdouble* get_acc_i();
        cdouble* get_acc_k();
        cdouble* get_dip_i();
        cdouble* get_dip_k();
        cdouble* get_pop();
        void write_diagnostics();
};

class DiagnosticsXZ: Diagnostics{
    protected:
        void _acc_i_buf(const int idx);
        void _acc_k_buf(const int idx);
        void _dip_i_buf(const int idx);
        void _dip_k_buf(const int idx);
        void _pop_buf(const int idx);
    public:
        DiagnosticsXZ();
        DiagnosticsXZ(Parameters* param);

};
#endif
