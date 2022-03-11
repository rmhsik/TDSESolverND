import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import h5py 
from mpi4py import MPI
from fileoperations import *
from pytdsesolver import param
from pytdsesolver import tdsesolver
import copy
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

filepath = 'data.txt'

def gendata(filepath):
    Es_arr = np.linspace(-0.1,0.1,10)
    phi_arr = np.linspace(0,2*np.pi,10)
    tmp = np.array(np.meshgrid(phi_arr,Es_arr)).T.reshape(-1,2)
    np.savetxt(filepath, tmp)

def main():
    data = np.empty(2,dtype=np.double)
    _param = copy.deepcopy(param.param)

    restart = False
    if len(sys.argv) > 1 and sys.argv[1] == 'restart': restart = True


    if rank == 0 and not restart:
        gendata(filepath)

    comm.barrier()
    while(getlinesnumber(filepath) > 0):
        if rank == 0:
            datafromfile = getdata(filepath,size)

            data = datafromfile[0][1]

            for i,val in enumerate(datafromfile[1:]):
                comm.Send(val[1], dest=i+1, tag=11)


        if rank != 0:
            comm.Recv(data,source=0,tag=11)

        comm.barrier()
        if (not np.isnan(data[0])) and (not np.isnan(data[1])):
            Esidx = str(format(data[1],".8f")).replace('.','')
            phiidx = str(format(data[0],".8f")).replace('.','')
            _param["phiBk"] = data[0]
            _param["E0j"] = data[1]
            _param["n_probes"] = 3
            _param["probe_def"] = (f"acc_i,results/acc_i_Es{Esidx}_Phi{phiidx}.dat;"
                                  f"acc_j,results/acc_j_Es{Esidx}_Phi{phiidx}.dat;"
                                  f"acc_k,results/acc_k_Es{Esidx}_Phi{phiidx}.dat")

            parameters = tdsesolver.Parameters()
            parameters.set(_param)
            parameters.print()
            tdse = tdsesolver.TDSESolver(parameters._obj)

            tdse.ipropagate()
            tdse.propagate()

            f = h5py.File(f'data/hhg_circular_Es{Esidx}_Phi{phiidx}.hdf5','a')
            acc_i = LoadComplexData(f'results/acc_i_Es{Esidx}_Phi{phiidx}.dat')
            acc_j = LoadComplexData(f'results/acc_j_Es{Esidx}_Phi{phiidx}.dat')
            acc_k = LoadComplexData(f'results/acc_k_Es{Esidx}_Phi{phiidx}.dat')
            f.create_dataset(f'acc_i', data=acc_i)
            f.create_dataset(f'acc_j', data=acc_j)
            f.create_dataset(f'acc_k', data=acc_k)
            f.create_dataset(f'phi', data=data[0])
            f.create_dataset(f'es', data=data[1])
            f.close()
            
        comm.barrier()
        if rank == 0:
            for i in range(size):
                removedataline(filepath,0)
        comm.barrier()

if __name__ == "__main__":
    main()

