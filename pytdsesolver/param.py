import numpy as np


w0 = 0.057
pi = np.pi
param = {"n_threads": 6,
         "init_wf": 1,
         "use_potential": 0,
         "geometry": 3,
         "ni": 600,
         "imin": -60.0,
         "imax": 60.0,
         "nj": 600,
         "jmin": -60.0,
         "jmax": 60.0,
         "nk": 600,
         "kmin": -60.0,
         "kmax": 60.0,
         "w0": w0,
         "period": 2*pi/w0,
         "tmax_ev": 4*2*pi/w0,
         "tmax_sim": 5*2*pi/w0,
         "dt": 0.02,
         "dt_ITP": 0.004,
         "nt": int(5*2*pi/w0/0.02),
         "nt_ITP": 2000,
         "nt_diag": 10,
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
         "B0k": 0.12,
         "phiEi": 0.0*pi,
         "phiEj": 0.0*pi,
         "phiEk": 0.0*pi,
         "phiBi": 0.0*pi,
         "phiBj": 0.0*pi,
         "phiBk": 0.0*pi,
         "n_probes": 2,
         "probe_def": ("acc_i,results/acc_i.dat;"
                       "acc_k,results/acc_k.dat")
        }
