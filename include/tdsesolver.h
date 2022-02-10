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

    class TDSESolver{
        private:
            Parameters *_param;

            Field *Afield_i;
            Field *Afield_j;
            Field *Afield_k;
            Field *Bfield_i;
            Field *Bfield_j;
            Field *Bfield_k;
            
            WF *_wf;
            Hamiltonian *_ham;
            Diagnostics *_diag;

            cdouble *_imask;
            cdouble *_kmask;

            double *_t, _dt;
            double *_i, _di;
            double *_j, _dj;
            double *_k, _dk;

            void (TDSESolver::*_propagate)();
            void (TDSESolver::*_ipropagate)();
            void _propagate_X();
            void _ipropagate_X();
            void _propagate_XZ();
            void _ipropagate_XZ();
            void _propagate_RZ();
            void _ipropagate_RZ();
            void _propagate_XYZ();
            void _ipropagate_XYZ();

            void _geom_X();
            void _geom_XZ();
            void _geom_RZ();
            void _geom_XYZ();
            void _fields_X();
            void _fields_XZ();
            void _fields_RZ();
            void _fields_XYZ();
            void _masks_X();
            void _masks_XZ();
            void _masks_RZ();
            void _masks_XYZ();
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
            void propagate();
            void ipropagate();
            
            ~TDSESolver();
    };

#endif
