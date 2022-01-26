#include <iostream>
#include <omp.h>
#include "tdsesolver.h"
#include "parameters.h"

int main(){
    Parameters param;
    param.n_threads = 4;

    param.geometry = XZ;
    param.init_wf = GAUS;
    param.ni = 100;
    param.nk = 100;
    param.imin = -100;
    param.imax =  100;
    param.kmin = -100;
    param.kmax =  100;
				
    param.w0 = 0.057;
    param.period = 2.0*M_PI/param.w0;
    param.tmax_ev = 4.0*param.period;
    param.tmax_sim = 5.0*param.period;
    param.dt = 0.02;
    param.dt_ITP = 0.004;
    param.nt = param.tmax_sim/param.dt;
    param.nt_ITP = 2000;
    param.nt_diag = 30;
			    
    param.env = SIN2;
    param.w0Ei = 0.057;
    param.w0Ek = 0.057;
    param.w0Bi = 0.057;
    param.w0Bk = 0.057;
								
    param.E0i = 0.067;
    param.E0k = 0.067;
    param.B0i = 0.0;
    param.B0k = 0.12;
													
    param.phiEi = 0.5*M_PI;
    param.phiEk = 0.0;
    param.phiBi = 0.0;
    param.phiBk = 0.0;

    char *acc_i_path = "results/acc_i.dat";
    param.set_acc_i_path(acc_i_path);
    char *acc_k_path = "results/acc_k.dat";
    param.set_acc_k_path(acc_k_path);
    char *dip_path = "results/dip8.dat";
    param.set_dip_path(dip_path); 
    //param.acc_path = "results/acc7.dat";
    //param.dip_path = "results/dip7.dat"; 

    param.print();
    TDSESolver tdse(&param);
    tdse.ipropagate();
    //tdse.propagate();
}
