#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cmath>
#include <string>

#define X    0
#define XZ   1
#define RZ   2

#define SIN2 0
#define TRAP 1

#define GAUS 0
#define EXPO 1

#define C  137.04

    class Parameters{
        private:
        public:
	    //Inital WF
	    int n_threads;
	    int init_wf;	 
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
            double w0Bi;     
            double w0Bk;     
            double E0i;      
            double E0k;      
            double B0i;      
            double B0k;      
            double phiEi;    
            double phiEk;    
            double phiBi;    
            double phiBk;
            int population;
            double pop_imin;
            double pop_imax;
            double pop_kmin;
            double pop_kmax;
            int acc_i;
            int acc_k;
	        std::string acc_i_path;
            std::string acc_k_path; 
	        std::string dip_i_path;
            std::string dip_k_path; 
            std::string pop_path;
        Parameters();
       	void check_param();	
	void print();
	void set_acc_i_path(char *val);
    void set_acc_k_path(char *val);
	void set_dip_i_path(char* val);
    void set_dip_k_path(char* val);
	void set_pop_path(char* val);
    };

#endif
