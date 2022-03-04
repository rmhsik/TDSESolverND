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


f = h5py.File('data/susceptibility_data_2.hdf5','a')

#Range of values:
B_array_scan = np.array([0.0,0.05,0.10,0.15,0.20])
E_array_scan = np.array([0.0340,0.0370, 0.0400, 0.0430, 0.0470, 0.0500])

#f = h5py.File(f'results/hhgcircular_phase_{idx}.hdf5','a')
#Range of values:
_param = param.param

for i in E_array_scan:
    count = 0
    for j in B_array_scan:
        _param["E0k"] = i
        Eidx = str(format(i,".5f")).replace('.','')
        Bidx = str(format(j,".5f")).replace('.','')
        _param["B0k"] = j 
        _param["probe_def"] = (f"dip_k,0.0,150.0,-150.0,150.0,results/dip_k_total_E{Eidx}_B{Bidx}.dat;"
                               f"dip_k,0.0,10.0,-10.0,10.0,results/dip_k_near_E{Eidx}_B{Bidx}.dat")

        parameters = tdsesolver.Parameters()
        parameters.set(_param)
        parameters.print()
        tdse = tdsesolver.TDSESolver(parameters._obj)

        tdse.ipropagate()
        tdse.propagate()
        tdse.delete()

        dip_k_total = LoadComplexData(f'results/dip_k_total_E{Eidx}_B{Bidx}.dat')
        dip_k_near  = LoadComplexData(f'results/dip_k_near_E{Eidx}_B{Bidx}.dat')
        Efield = np.genfromtxt('results/Efield_k.dat')

        dataset_Eprefix = 'E_'+str(format(i,".5f")).replace('.','')
        dataset_Bprefix = 'B_'+str(format(j,".5f")).replace('.','')
        f.create_dataset(f'{dataset_Eprefix}/{dataset_Bprefix}/dip_k_total',data=dip_k_total)
        f.create_dataset(f'{dataset_Eprefix}/{dataset_Bprefix}/dip_k_near',data=dip_k_near)
        if count==0:
            f.create_dataset(f'{dataset_Eprefix}/Efield',data=Efield)
            count=1
      

f.close()
