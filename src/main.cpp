#include <iostream>
#include <omp.h>
#include <vector>
#include "tdsesolver.h"
#include "parameters.h"
#include "diagnostics.h"

int main(){
    Parameters param;
    param.n_threads = 4;

    param.geometry = XYZ;
    param.init_wf = EXPO;
    param.use_potential = 0;
    param.ni = 10;
    param.nk = 10;
    param.nj = 10;
    param.imin = -20;
    param.imax = 20;
    param.jmin = 20;
    param.jmax = 20;
    param.kmin = -20;
    param.kmax =  20;
				
    param.w0 = 0.057;
    param.period = 2.0*M_PI/param.w0;
    param.tmax_ev = 4.0*param.period;
    param.tmax_sim = 5.0*param.period;
    param.dt = 0.04;
    param.dt_ITP = 0.004;
    param.nt = 20;//param.tmax_sim/param.dt;
    param.nt_ITP = 20;
    param.nt_diag = 1;
			    
    param.env = SIN2;
    param.w0Ei = 0.057;
    param.w0Ek = 0.057;
    param.w0Bi = 0.057;
    param.w0Bk = 0.057;
								
    param.E0i = 0.067;
    param.E0j = 0.010;
    param.E0k = 0.067;
    param.B0i = 0.000;
    param.B0k = 0.12;
													
    param.phiEi = 0.5*M_PI;
    param.phiEk = 0.0;
    param.phiBi = 0.0;
    param.phiBk = 0.00;

    param.n_probes = 3;
    param.probe_def = "acc_i,results/acc_i.dat;acc_j,results/acc_j.dat;acc_k,results/acc_k.dat";
    param.check_param();
    param.print();
    TDSESolver *tdse;
    tdse = new TDSESolver(&param);
    tdse->ipropagate();
    tdse->propagate();
    delete tdse;
}
