#include <iostream>
#include <omp.h>
#include <vector>
#include "tdsesolver.h"
#include "parameters.h"
#include "diagnostics.h"

int main(){
    Parameters param;
    param.n_threads = 6;

    param.geometry = RZ;
    param.init_wf = GAUS;
    param.use_potential = 0;
    param.ni = 600;
    param.nk = 1000;
    param.imin =  -120;
    param.imax =  120;
    param.kmin = -120;
    param.kmax =  120;
				
    param.w0 = 0.057;
    param.period = 2.0*M_PI/param.w0;
    param.tmax_ev = 4.0*param.period;
    param.tmax_sim = 5.0*param.period;
    param.dt = 0.02;
    param.dt_ITP = 0.004;
    param.nt = param.tmax_sim/param.dt;
    param.nt_ITP = 2000;
    param.nt_diag = 100;
			    
    param.env = SIN2;
    param.w0Ei = 0.057;
    param.w0Ek = 0.057;
    param.w0Bi = 0.057;
    param.w0Bk = 0.057;
								
    param.E0i = 0.000;
    param.E0k = 0.067;
    param.B0i = 0.000;
    param.B0k = 0.12;
													
    param.phiEi = 0.0*M_PI;
    param.phiEk = 0.0;
    param.phiBi = 0.0;
    param.phiBk = 0.00;

    param.n_probes = 2;
    param.probe_def = "acc_i,results/acc_i.dat;acc_k,results/acc_k.dat";
    param.check_param();
    param.print();
    TDSESolver *tdse;
    tdse = new TDSESolver(&param);
    tdse->ipropagate();
    tdse->propagate();
    delete tdse;
}
