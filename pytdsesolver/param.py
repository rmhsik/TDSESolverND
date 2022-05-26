import numpy as np


w0 = 0.057
pi = np.pi
param = {"n_threads": 16,
         "init_wf": 2,
         "use_potential": 2,
         "geometry": 3,
         "ni": 800,
         "imin": -100.0,
         "imax": 10.0,
         "nj": 800,
         "jmin": -100.0,
         "jmax": 100.0,
         "nk": 800,
         "kmin": -100.0,
         "kmax": 100.0,
         "w0": w0,
         "period": 2*pi/w0,
         "tmax_ev": 4*2*pi/w0,
         "tmax_sim": 4*2*pi/w0,
         "dt": 0.02,
         "dt_ITP": 0.01,
         "nt": int(4*2*pi/w0/0.02),
         "nt_ITP": 2000,
         "nt_diag": 100,
         "env": 0,
         "w0Ei": w0,
         "w0Ej": w0, 
         "w0Ek": w0,
         "w0Bi": w0,
         "w0Bj": w0,
         "w0Bk": w0,
         "E0i": 0.000,
         "E0j": 0.000,
         "E0k": 0.067,
         "B0i": 0.00,
         "B0j": 0.00,
         "B0k": 0.00,
         "phiEi": 0.5*pi,
         "phiEj": 0.0*pi,
         "phiEk": 0.0*pi,
         "phiBi": 0.0*pi,
         "phiBj": 0.0*pi,
         "phiBk": 0.0*pi,
         "n_probes": 3,
         "probe_def": ("acc_i,results/acc_i_1.dat;"
                       "acc_j,results/acc_j_1.dat;"
                       "acc_k,results/acc_k_1.dat;"
                       )
        }
