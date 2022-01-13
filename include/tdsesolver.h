#ifndef TDSESOLVER_H
#define TDSESOLVER_H

#define cdouble std::complex<double>
#define I std::complex<double>(0.0,1.0)

#include <complex>
#include "parameters.h"
#include "fields.h"
#include "wavefunction.h"
#include "hamiltonian.h"

    class TDSESolver{
        private:
            Parameters _param;

            Field Afield_i;
            Field Afield_k;
            Field Bfield_i;
            Field Bfield_k;
            
            WF _wf;
            Hamiltonian _ham;

            cdouble *_accmask;
            cdouble *_imask;
            cdouble *_kmask;

            double *_t, _dt;
            double *_i, _di;
            double *_k, _dk;

	    void (TDSESolver::*_propagate)();
	    void (TDSESolver::*_ipropagate)();

        public:
            TDSESolver();
            TDSESolver(Parameters param);
            void setup_time();
            void setup_geometry();
            void setup_fields();
            void setup_masks();
            void setup_wf();
            void setup_ham();
	    void propagate();
	    void ipropagate();
            void propagate_X();
            void ipropagate_X();
	    void propagate_RZ();
	    void ipropagate_RZ();
    };


#endif
