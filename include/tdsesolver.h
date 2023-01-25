#ifndef TDSESOLVER_H
#define TDSESOLVER_H

#define cdouble std::complex<double>
#define I std::complex<double>(0.0,1.0)

#include <complex>
#include "parameters.h"
#include "fields.h"
#include "wavefunction.h"
#include "hamiltonian.h"
#include "diagnostics.h"
#include "mpi_grid.h"

    class TDSESolver{
        private:
            Parameters *_param;

            mpi_grid *_mpi_grid;


            Field *Afield_i;
            Field *Afield_k;
            
            WF *_wf;
            Hamiltonian *_ham;
            Diagnostics *_diag;

            cdouble *_imask;
            cdouble *_kmask;

            double *_t, _dt;
            double *_i, _di;
            double *_k, _dk;

            void (TDSESolver::*_propagate)();
            void (TDSESolver::*_ipropagate)();
            void _propagate_XZ();
            void _ipropagate_XZ();

            void _geom_XZ();
            void _fields_XZ();
            void _masks_XZ();

        public:
            TDSESolver();
            TDSESolver(Parameters *param);
            void setup_time();
            void setup_geometry();
            void setup_fields();
            void setup_masks();
            void setup_diagnostics();
            void setup_wf();
            void setup_ham();
            void setup_mpi();
            void propagate();
            void ipropagate();
            
            ~TDSESolver();
    };

#endif
