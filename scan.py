import numpy as np
from pytdsesolver import param
from pytdsesolver import tdsesolver

#Range of values:
B_array_scan = np.array([0.0,0.05,0.10,0.15,0.20])
B_array_scan = np.array([0.0])
E0 = 0.0400

_param = param.param
parameters = tdsesolver.Parameters()
for i in B_array_scan:
    _param["E0k"] = E0
    _param["B0k"] = i
    E0_name = "{:.4f}".format(E0).replace('.','')
    B0_name = "{:.4f}".format(i).replace('.','')
    _param["acc_path"] = f"results/acc_E{E0_name}_B{B0_name}.dat"
    _param["dip_path"] = f"results/dip_E{E0_name}_B{B0_name}.dat"
    _param["pop_path"] = f"results/pop_E{E0_name}_B{B0_name}.dat"
    parameters.set(_param)
    parameters.print()
    tdse = tdsesolver.TDSESolver(parameters._obj)
    tdse.ipropagate()
    tdse.propagate()
