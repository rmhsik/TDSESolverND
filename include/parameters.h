#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cmath>
#include <string>

#define XZ   1

#define SIN2 0
#define TRAP 1

#define GAUS 0
#define EXPO 1
#define GAUS_ANTI 2

#define USER_DEFINED 0
#define HYDROGEN 1
#define ARGON 2
#define C  137.04

    class Parameters{
        private:
        public:
	    //Inital WF
	    int n_threads;
	    int init_wf;	 
            int use_potential;
            int geometry;           
            int ni;  
            double imin;
            double imax;     
            int nk;         
            double kmin;     
            double kmax;     
            double w0;       
            double period;   
            double tmax_ev;  
            double tmax_sim; 
            double dt;       
            double dt_ITP;   
            int nt;          
            int nt_ITP;      
            int nt_diag ;    
            int env;         
            double w0Ei;     
            double w0Ek;    
            double E0i;      
            double E0k;
            double phiEi;    
            double phiEk;
            int n_probes;
            std::string probe_def;
            Parameters();

            void check_param();	
            void print();
            void set_probe_def(char *val);
    };

#endif
