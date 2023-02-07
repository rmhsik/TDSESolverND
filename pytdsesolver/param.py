import numpy as np


w0 = 0.057/3
pi = np.pi
dt = 0.04
param = {"n_threads": 16,
         "init_wf": 1, #EXPO
         "use_potential": 1,
         "geometry": 1, #XZ
         "ni": 2048,
         "imin": -256.0,
         "imax": 256.0,
         "nk": 4096,
         "kmin": -512.0,
         "kmax": 512.0,
         "w0": w0,
         "period": 2*pi/w0,
         "tmax_ev": 4*2*pi/w0,
         "tmax_sim": 5*2*pi/w0,
         "dt": dt,
         "dt_ITP": 0.005,
         "nt": int(5*2*pi/w0/dt),
         "nt_ITP": 5000,
         "nt_diag": 100,
         "env": 0, #SIN2
         "w0Ei": w0,
         "w0Ek": w0,
         "E0i": 0.000,
         "E0k": 0.067,
         "phiEi": 0.0*pi,
         "phiEk": 0.0*pi,
         "kk": w0/137.04,
         "n_probes": 2,
         "probe_def": ("acc_i,results/acc_i_3_ND.dat;"
                       "acc_k,results/acc_k_3_ND.dat;"
                       )
        }
