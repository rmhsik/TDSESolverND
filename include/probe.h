#ifndef PROBE_H
#define PROBE_H
#include "wavefunction.h"
#include "hamiltonian.h"
#include "parameters.h"


class Probe{
    protected:
        cdouble *_data;
        std::string _def;
        
        WF *_wf;
        Hamiltonian *_ham;
        Parameters *_param;
        cdouble *_tempmask;
        double *_i, *_j, *_k;
        int _ni, _nj, _nk, _nt, _nt_diag;
        double _di, _dj, _dk;

        int _type;
        std::string _data_path;
        double _int_imin, _int_imax;
        double _int_kmin, _int_kmax;
        double _int_jmin, _int_jmax;
        void _parse_def();

        virtual void _acc_i(const int idx){};
        virtual void _acc_k(const int idx){};
        virtual void _dip_i(const int idx){};
        virtual void _dip_k(const int idx){};
        virtual void _pop(const int idx){};
        virtual void _dens(const int idx){};
        virtual void _wf_snap(const int idx){};
        void (Probe::*_calc)(const int idx);
    public:
        Probe();
        Probe(std::string def);

        void set_wf(WF *wf);
        void set_ham(Hamiltonian *ham);
        void set_param(Parameters *param);
        void set_geometry(double* i, double *j, double* k, const double di, const double dj, const double dk);
        void set_tempmask(cdouble* tempmask);

        void write_probe();
        cdouble* get_data();
        void calc(const int idx);
        
        ~Probe();
};


#endif
