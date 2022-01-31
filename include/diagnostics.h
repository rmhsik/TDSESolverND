#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include <string>
#include <vector>
#include "wavefunction.h"
#include "hamiltonian.h"
#include "parameters.h"
#include "probe_factory.h"

#define ACC_I 0
#define ACC_K 1
#define DIP_I 2
#define DIP_K 3
#define POP   4

class Diagnostics{
    private:
        std::vector<Probe*> _probe_vec;
        std::vector<std::string> _def_vec;
        int _n_probes;

        Parameters *_param;
        cdouble *_tempmask;
        
        int _ni, _nk;
        double _di, _dk; 
        double *_i, *_k, *_t;
        Hamiltonian *_ham;
        WF *_wf;
        

                
    public:
        Diagnostics();
        Diagnostics(int n_probes, std::string def);
        void run_diagnostics(const int idx);
        void write_diagnostics();
        void set_parameters(Parameters *param);
        void set_geometry(double *i, double *k, double *t, const double di, const double dk);
        void set_ham(Hamiltonian *ham);
        void set_wf(WF *wf); 
        void set_tempmask();
        void create_probes();

};

#endif
