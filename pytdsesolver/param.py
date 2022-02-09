import numpy as np


w0 = 0.057
pi = np.pi
param = {"n_threads": 6,
         "init_wf": 1,
         "use_potential": 1,
         "geometry": 0,
         "ni": 10000,
         "imin": -200.0,
         "imax": 200.0,
         "nk": 500,
         "kmin": -120.0,
         "kmax": 120.0,
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
         "w0Ek": w0,
         "w0Bi": w0,
         "w0Bk": w0,
         "E0i": 0.000,
         "E0k": 0.067,
         "B0i": 0.00,
         "B0k": 0.12,
         "phiEi": 0.0*pi,
         "phiEk": 0.5*pi,
         "phiBi": 0.0,
         "phiBk": -0.0*pi,
         "n_probes": 2,
         "probe_def": ("dens,results/dens5.dat;"
                       "wf,results/wf_snap5.dat")
        }
