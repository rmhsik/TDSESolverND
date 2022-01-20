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

const int n_threads = 6;
    class Parameters{
        private:
        public:
	    //Inital WF
	    int init_wf	    = GAUS;
            //Geometry
            int geometry    = XZ;            
            int ni          = 1000;
            double imin     = -100.0;
            double imax     = 100.0;
            int nk          = 600;
            double kmin     = -120.0;
            double kmax     = 120.0;

            // Time
            double w0       = 0.057;
            double period   = 2.0*M_PI/w0;
            double tmax_ev  = 4.0*period;
            double tmax_sim = 5.0*period;
            double dt       = 0.02;
            double dt_ITP   = 0.004;
            int nt          = tmax_sim/dt;
            int nt_ITP      = 2000;
            int nt_diag     = 100;

            //Fields
            int env         = SIN2;
            double w0Ei     = w0;
            double w0Ek     = w0;
            double w0Bi     = w0;
            double w0Bk     = w0;
            
            double E0i      = 0.067;
            double E0k      = 0.067;
            double B0i      = 0.0;
            double B0k      = 0.12;

            double phiEi    = 0.5*M_PI;
            double phiEk    = 0.0;
            double phiBi    = 0.0;
            double phiBk    = 0.0;

	    //File paths
	    std::string acc_path = "results/acc2.dat";
	    std::string dip_path = "results/dip2.dat";
        Parameters(){
            if(geometry==X){
                nk = 1;
                kmin = 0.0;
                kmax = 0.0;
            }
            else if(geometry == RZ){
                imin = 0.0;
            }
        }       
    };

#endif
