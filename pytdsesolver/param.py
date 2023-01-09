import numpy as np


w0 = 0.057/3
pi = np.pi
param = {"n_threads": 16,
         "init_wf": 1, #EXPO
         "use_potential": 1,
         "geometry": 1, #XZ
         "ni": 3072,
         "imin": -2000.0,
         "imax": 2000.0,
         "nk": 6144,
         "kmin": -2000.0,
         "kmax": 2000.0,
         "w0": w0,
         "period": 2*pi/w0,
         "tmax_ev": 4*2*pi/w0,
         "tmax_sim": 5*2*pi/w0,
         "dt": 0.02,
         "dt_ITP": 0.005,
         "nt": int(5*2*pi/w0/0.02),
         "nt_ITP": 5000,
         "nt_diag": 100,
         "env": 0, #SIN2
         "w0Ei": w0,
         "w0Ek": w0,
         "E0i": 0.000,
         "E0k": 0.067,
         "phiEi": 0.0*pi,
         "phiEk": 0.0*pi,
         "n_probes": 2,
         "probe_def": ("acc_i,results/acc_i_3.dat;"
                       "acc_k,results/acc_k_3.dat;"
                       )
        }
