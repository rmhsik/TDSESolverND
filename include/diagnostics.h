#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include <string>
#include <vector>
#include "wavefunction.h"
#include "hamiltonian.h"
#include "parameters.h"
#include "probe_factory.h"
#include "mpi_grid.h"

class Diagnostics{
    private:
        std::vector<Probe*> _probe_vec;
        std::vector<std::string> _def_vec;
        int _n_probes;

        mpi_grid *_mpi_grid;
        Parameters *_param;
        cdouble *_tempmask;
        
        int _ni,_nk;
        double _di,_dk; 
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
        void set_mpi(mpi_grid *grid);
        void create_probes();

        ~Diagnostics();

};

#endif
