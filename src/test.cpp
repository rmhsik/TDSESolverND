#include <iostream>
#include <cmath>
#include <complex>
#include <mpi.h>

int main(){
    MPI_Init(NULL,NULL);

    int size_process;
    MPI_Comm_size(MPI_COMM_WORLD,&size_process);
    int num_process;
    MPI_Comm_rank(MPI_COMM_WORLD,&num_process);
    std::cout<<"World size: "<<size_process<<" Num: "<<num_process<<std::endl;

    return 0;
}
