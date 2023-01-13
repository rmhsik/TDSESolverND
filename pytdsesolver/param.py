import numpy as np


w0 = 0.057/3.0
pi = np.pi
param = {"n_threads": 16,
         "init_wf": 1, #EXPO
         "use_potential": 1,
         "geometry": 1, #XZ
         "ni": 512,
         "imin": -256.0,
         "imax": 256.0,
         "nk": 16384,
         "kmin": -512.0,
         "kmax": 512.0,
         "w0": w0,
         "period": 2*pi/w0,
         "tmax_ev": 4*2*pi/w0,
         "tmax_sim": 5*2*pi/w0,
         "dt": 0.08,
         "dt_ITP": 0.005,
         "nt": int(5*2*pi/w0/0.08),
         "nt_ITP": 5000,
         "nt_diag": 100,
         "env": 0, #SIN2
         "w0Ei": w0,
         "w0Ek": w0,
         "E0i": 0.000,
         "E0k": 0.067,
         "phiEi": 0.0*pi,
         "phiEk": 0.0*pi,
         "kk": 0.0001386456509048453,
         "n_probes": 2,
         "probe_def": ("acc_i,results/acc_i_3_ND_E0.067.dat;"
                       "acc_k,results/acc_k_3_ND_E0.067.dat;"
                       )
        }
