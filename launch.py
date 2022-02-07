from pytdsesolver import param
from pytdsesolver import tdsesolver

_param = param.param
parameters = tdsesolver.Parameters()
parameters.set(_param)
parameters.print()
tdse = tdsesolver.TDSESolver(parameters._obj)
#tdse.ipropagate()
tdse.propagate()
