#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include <string>
#include "wavefunction.h"
#include "hamiltonian.h"
#include "parameters.h"

#define ACC_I 0
#define ACC_K 1
#define DIP_I 2
#define DIP_K 3
#define POP   4

class Diagnostics{
    protected:
        Parameters *_param;
        cdouble *_data;
        
        std::sttring _def;
        int _type;
        std::string _data_path;
        double _int_imin, _int_imax;
        double _int_kmin, _int kmax;

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

        void _parse_def();
    public:
        Diagnostics();
        Diagnostics(std::string def);
        void set_parameters(Parameters *param);
        void set_geometry(double *i, double *k, const double di, const double dk);
        void set_ham(Hamiltonian *ham);
        void set_wf(WF *wf); 
        void run_diagnostics(const int idx);
        cdouble* get_data();
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
