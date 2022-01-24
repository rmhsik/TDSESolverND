import param
import tdsesolver

_param = param.param
parameters = tdsesolver.Parameters()
parameters.set(_param)
parameters.print()
tdse = tdsesolver.TDSESolver(parameters._obj)
tdse.ipropagate()
