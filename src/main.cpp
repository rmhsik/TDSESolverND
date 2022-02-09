#include <iostream>
#include <omp.h>
#include <vector>
#include "tdsesolver.h"
#include "parameters.h"
#include "diagnostics.h"

int main(){
    Parameters param;
    param.n_threads = 6;

    param.geometry = X;
    param.init_wf = GAUS;
    param.use_potential = 0;
    param.ni = 5000;
    param.nk = 1;
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
								
    param.E0i = 0.067;
    param.E0k = 0.067;
    param.B0i = 0.000;
    param.B0k = 0.012;
													
    param.phiEi = 0.5*M_PI;
    param.phiEk = 0.0;
    param.phiBi = 0.0;
    param.phiBk = 0.0;

    param.n_probes = 1;
    param.probe_def = "dens,results/dens.dat;";

    //param.acc_path = "results/acc7.dat";
    //param.dip_path= "results/dip7.dat"; 
    param.check_param();
    param.print();
    for(int i=0;i<20;i++){
        TDSESolver *tdse;
        tdse = new TDSESolver(&param);
        //tdse.ipropagate();
        //tdse.propagate();
        delete tdse;
    }
}
