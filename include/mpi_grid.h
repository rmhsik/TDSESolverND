#ifndef MPI_GRID_H
#define MPI_GRID_H

#include <mpi.h>


struct mpi_grid{
    MPI_Comm comm;
    int dims[2];
    int coords[2];
    int rank;
    int size;
    int reorder = 1;
    int period [2];

};

#endif
