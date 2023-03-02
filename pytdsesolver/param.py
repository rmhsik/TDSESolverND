import numpy as np

w0 = 0.0117
pi = np.pi
c = 137.04
dt = 0.16

param = {"n_threads": 72,
         "init_wf": 1, #EXPO
         "use_potential": 1,
         "geometry": 1, #XZ
         "ni": 2**12,
         "imin": -1024.0,
         "imax": 1024.0,
         "nk": 2**16,
         "kmin": -4096.0,
         "kmax": 4096.0,
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
         "probe_def": ("acc_i,results/acc_i_3900_D_E0.067_2.dat;"
                       "acc_k,results/acc_k_3900_D_E0.067_2.dat;"
                       )
        }
