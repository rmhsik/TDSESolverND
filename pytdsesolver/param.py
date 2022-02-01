param = {"n_threads": 6,
         "init_wf": 1,
         "geometry": 2,
         "ni": 500,
         "imin": -100.0,
         "imax": 100.0,
         "nk": 1400,
         "kmin": -120.0,
         "kmax": 120.0,
         "w0": 0.057,
         "period": 2*3.141592/0.057,
         "tmax_ev": 4*2*3.141592/0.057,
         "tmax_sim": 5*2*3.141592/0.057,
         "dt": 0.02,
         "dt_ITP": 0.004,
         "nt": int(5*2*3.141592/0.057/0.02),
         "nt_ITP": 2000,
         "nt_diag": 50,
         "env": 0,
         "w0Ei": 0.057,
         "w0Ek": 0.057,
         "w0Bi": 0.057,
         "w0Bk": 0.057,
         "E0i": 0.000,
         "E0k": 0.067,
         "B0i": 0.00,
         "B0k": 0.12,
         "phiEi": 0.0*3.141592,
         "phiEk": 0.0,
         "phiBi": 0.0,
         "phiBk": 0.0,
         "n_probes": 2,
         "probe_def": ("acc_i,results/acc_i.dat;"
                        "acc_k,results/acc_k.dat")
        }
