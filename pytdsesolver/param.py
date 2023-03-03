import numpy as np

w0 = 0.057
pi = np.pi
dt = 0.02
param = {"n_threads": 16,
         "init_wf": 1, #EXPO
         "use_potential": 2,
         "geometry": 1, #XZ
         "ni": 2048,
         "imin": -512.0,
         "imax": 512.0,
         "nk": 4096,
         "kmin": -1200.0,
         "kmax": 1200.0,
         "w0": w0,
         "period": 2*pi/w0,
         "tmax_ev": 4*2*pi/w0,
         "tmax_sim": 5*2*pi/w0,
         "dt": dt,
         "dt_ITP": 0.005,
         "nt": int(5*2*pi/w0/dt),
         "nt_ITP": 100,
         "nt_diag": 10,
         "env": 0, #SIN2
         "w0Ei": w0,
         "w0Ek": w0,
         "E0i": 0.000,
         "E0k": 0.067,
         "phiEi": 0.0*pi,
         "phiEk": 0.0*pi,
         "n_probes": 2,
         "probe_def": ("acc_i,results/acc_i_800_ND_helium.dat;"
                       "acc_k,results/acc_k_800_ND_helium.dat;"
                       )
        }
