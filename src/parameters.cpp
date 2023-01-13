#include <iostream>
#include <cmath>
#include <string>
#include "parameters.h"

Parameters::Parameters(){
        n_threads = 6;
        //Inital WF
        init_wf	    = GAUS;

        // Use potential
        use_potential = 0;

        //Geometry
        geometry    = XZ;            
        ni          = 500;
        imin        = -100.0;
        imax        = 100.0;
        nk          = 500;
        kmin        = -120.0;
        kmax        = 120.0;

        // Time
        w0          = 0.057;
        period      = 2.0*M_PI/w0;
        tmax_ev     = 4.0*period;
        tmax_sim    = 5.0*period;
        dt          = 0.02;
        dt_ITP      = 0.004;
        nt          = tmax_sim/dt;
        nt_ITP      = 2000;
        nt_diag     = 100;

        //Fields
        env         = SIN2;
        w0Ei        = w0;
        w0Ek        = w0;
        
        E0i         = 0.067;
        E0k         = 0.067;

        phiEi       = 0.5*M_PI;
        phiEk       = 0.0;
        
        kk          = 0.0; 
        // Diagnostics selection and parameters
        n_probes    = 2;

	    //File paths
        probe_def   = "acc_i,results/acc_i_test_diag.dat;"
                      "acc_k,results/acc_k_test_diag.dat";
        check_param();
}

void Parameters::check_param(){
}

void Parameters::set_probe_def(char *val){
    probe_def = std::string(val);
}


void Parameters::print(){
    std::cout<<"Parameters:\n";
    std::cout<<"------------------------\n";
    std::cout<<"\tn_threads: "<<n_threads<<std::endl;
    std::cout<<"\tinit_wf: "<<init_wf<<std::endl;
    std::cout<<"\tuse_potential: "<<use_potential<<std::endl;
    std::cout<<"\tgeometry: "<<geometry<<std::endl;
    std::cout<<"\tni: "<<ni<<std::endl;
    std::cout<<"\timin: "<<imin<<std::endl;
    std::cout<<"\timax: "<<imax<<std::endl;
    std::cout<<"\tnk: "<<nk<<std::endl;
    std::cout<<"\tkmin: "<<kmin<<std::endl;
    std::cout<<"\tkmax: "<<kmax<<std::endl;
    std::cout<<"\tw0: "<<w0<<std::endl;
    std::cout<<"\tperiod: "<<period<<std::endl;
    std::cout<<"\ttmax_ev: "<<tmax_ev<<std::endl;
    std::cout<<"\ttmax_sim: "<<tmax_sim<<std::endl;
    std::cout<<"\tdt: "<<dt<<std::endl;
    std::cout<<"\tdt_ITP: "<<dt_ITP<<std::endl;
    std::cout<<"\tnt: "<<nt<<std::endl;
    std::cout<<"\tnt_ITP: "<<nt_ITP<<std::endl;
    std::cout<<"\tnt_diag: "<<nt_diag<<std::endl;
    std::cout<<"\tenv: "<<env<<std::endl;
    std::cout<<"\tw0Ei: "<<w0Ei<<std::endl;
    std::cout<<"\tw0Ek: "<<w0Ek<<std::endl;
    std::cout<<"\tE0i: "<<E0i<<std::endl;
    std::cout<<"\tE0k: "<<E0k<<std::endl;
    std::cout<<"\tphiEi: "<<phiEi<<std::endl;
    std::cout<<"\tphiEk: "<<phiEk<<std::endl;
    std::cout<<"\tkk: "<<kk<<std::endl;
    std::cout<<"\tn_probes: "<<n_probes<<std::endl;
    std::cout<<"\tprobe_def: "<<probe_def<<std::endl;
    std::cout<<std::endl;
}
