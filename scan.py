import numpy as np
import h5py
from pytdsesolver import param
from pytdsesolver import tdsesolver

def LoadComplexData(file,**genfromtext_args):
        """
        Load complex data in the C++ format in numpy.
        """
        array_as_strings = np.genfromtxt(file,dtype=str,**genfromtext_args)
        complex_parser = np.vectorize(lambda x: complex(*eval(x)))
        return complex_parser(array_as_strings)

idx = 0
phi_min = 0.0
phi_max = 0.5

#f = h5py.File(f'results/hhgcircular_phase_{idx}.hdf5','a')
#Range of values:
phiBk = np.linspace(phi_min,phi_max, 100)
count = 13
_param = param.param

for i in phiBk:
    _param["phiBk"] = i*np.pi
    _param["probe_def"] = f"acc_i,results/acc_i_temp_{idx}.dat;acc_k,results/acc_k_temp_{idx}.dat"
    parameters = tdsesolver.Parameters()
    parameters.set(_param)
    parameters.print()
    tdse = tdsesolver.TDSESolver(parameters._obj)
    #tdse.ipropagate()
    #tdse.propagate()
    #del tdse
    tdse.delete()
    #del tdse
    count += 1
    #acc_i = LoadComplexData(f'results/acc_i_temp_{idx}.dat')
    #acc_k = LoadComplexData(f'results/acc_k_temp_{idx}.dat')
    #dataset_prefix = 'phi_'+str(format(i,".5f")).replace('.','')
    #f.create_dataset(f'{dataset_prefix}/acc_i',data=acc_i)
    #f.create_dataset(f'{dataset_prefix}/acc_k',data=acc_k)

#f.close()
