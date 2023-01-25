#include <iostream>
#include <tuple>
#include <string>
#include <cstdint>
#include <cstring>
#include <omp.h>
#include "tdsesolver.h"
#include "debug.h"
#include "utils.h"

TDSESolver::TDSESolver(){
}

TDSESolver::TDSESolver(Parameters *param){
    std::string path;
    _param = param;
    _param->check_param();

    setup_time();

    omp_set_num_threads(_param->n_threads);
    setup_mpi();
    MPI_Barrier(_mpi_grid->comm);
    if(_mpi_grid->rank==0){
        _param->print();
    }
    setup_geometry();
    setup_wf();
    setup_fields();
    setup_ham();
    setup_masks();
    setup_diagnostics();
    }

void TDSESolver::setup_mpi(){
    MPI_Init(NULL,NULL);
    _mpi_grid = new mpi_grid;
    _mpi_grid->dims[0] = 3;
    _mpi_grid->dims[1] = 2;
    _mpi_grid->period[0] = 0;
    _mpi_grid->period[1] = 0;
    _mpi_grid->reorder = 0;

    MPI_Comm_size(MPI_COMM_WORLD,
                  &_mpi_grid->size);
    MPI_Dims_create(_mpi_grid->size,2,_mpi_grid->dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2,
                    _mpi_grid->dims,
                    _mpi_grid->period,
                    _mpi_grid->reorder,
                    &_mpi_grid->comm);
    MPI_Comm_rank(_mpi_grid->comm,&_mpi_grid->rank);
    MPI_Cart_coords(_mpi_grid->comm, _mpi_grid->rank, 2, _mpi_grid->coords);
    if(_mpi_grid->rank==0){ 
        std::cout<<"[MPI] Size: "<<_mpi_grid->size<<std::endl;
    }
    std::cout<<"[MPI] Node number: "<<_mpi_grid->rank<<". Coords=("<<_mpi_grid->coords[0]<<","<<_mpi_grid->coords[1]<<")"<<std::endl;

}

void TDSESolver::setup_time(){
    std::string path;
    std::tie(_t,_dt) = linspace<double>(0.0,_param->tmax_sim, _param->nt);
    path = "results/time.dat";
    write_array(_t,_param->nt,path);
}

void TDSESolver::setup_geometry(){
    std::string path;
    switch(_param->geometry){
        case XZ:
            _geom_XZ();
    }
    path = "results/i.dat";
    write_array(_i,_param->ni,path);
    path = "results/k.dat";
    write_array(_k,_param->nk,path);
}

void TDSESolver::setup_fields(){
    switch(_param->geometry){
        case XZ:
            _fields_XZ();
            break;
    }
}

void TDSESolver::setup_wf(){
    _wf = new WF(_param);
    _wf->set_mpi(_mpi_grid);
    _wf->set_geometry(_i,_k,_di,_dk);
    switch(_param->init_wf){
    	case GAUS:
            _wf->gaussian(0.0,0.0,1.0);
            break;
        case EXPO:
            _wf->exponential(0.0,0.0,1.0);
            break;
    } 
    cdouble norm = _wf->norm();
    double tstart, tend;
    if(_mpi_grid->rank==0){
        std::cout<<"norm before: "<<norm<<std::endl;
    }
    (*_wf) /= norm;
    tstart = omp_get_wtime();
    norm = _wf->norm();
    tend = omp_get_wtime();
    if(_mpi_grid->rank==0){
        std::cout<<"norm after: "<<norm<<" Time to calc norm: "<<tend-tstart<<std::endl;
    }
}

void TDSESolver::setup_ham(){
    _ham = new Hamiltonian(_param);
    _ham->set_geometry(_i,_k,_t,_di,_dk,_dt);
    _ham->set_fields(Afield_k);
}

void TDSESolver::setup_masks(){
    _imask = new cdouble[_param->ni];
    _kmask = new cdouble[_param->nk];
    std::string path;    
    
    switch(_param->geometry){
        case XZ:
            _masks_XZ();
            break; 
    }

}

void TDSESolver::setup_diagnostics(){
    _diag = new Diagnostics(_param->n_probes, _param->probe_def);
    _diag->set_mpi(_mpi_grid);
    _diag->set_parameters(_param);
    _diag->set_geometry(_i,_k,_t,_di,_dk);
    _diag->set_ham(_ham);
    _diag->set_wf(_wf);
    _diag->set_tempmask();
    _diag->create_probes();
}

void TDSESolver::ipropagate(){
    (this->*(this->_ipropagate))();
	//_wf->save_wf2();
    //_wf->set_to_ground();
}

void TDSESolver::propagate(){
    (this->*(this->_propagate))();
}

TDSESolver::~TDSESolver(){
    MPI_Barrier(_mpi_grid->comm);
    std::cout<<"[MPI] Deleting node: "<<_mpi_grid->rank<<std::endl;
    delete[] _t;
    delete[] _i;
    delete[] _k;
    delete[] _imask;
    delete[] _kmask;
    delete Afield_i;
    delete Afield_k;
    delete _wf;
    delete _ham;
    delete _diag;
    MPI_Finalize();
}
