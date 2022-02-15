#include <iostream>
#include <cmath>
#include <complex>
#include <mpi.h>

double gauss(double i, double j, double k){
    return exp(-i*i-j*j-k*k);
}

int main(){
    MPI_Init(NULL,NULL);
    MPI_Status status;

    int size_process;
    MPI_Comm_size(MPI_COMM_WORLD,&size_process);
    
    double *i, *j, *k, ***wf;
    int nmax = 8;
    int ni = 1000;
    int nj = 500;
    int nk = 100;
    double imax= 50.0;
    double jmax = 50.0;
    double kmax = 50.0;
    double di = 2.0*imax/ni;
    double dk = 2.0*kmax/nk;
    double dj = 2.0*jmax/nj;
    double norm;
    i = new double[ni];
    j = new double[nj];
    k = new double[nk];

    wf = new double**[ni/size_process];
    for(int ii=0; ii<ni/size_process; ii++){
        wf[ii] = new double*[nj];
        for(int jj=0;jj<nj;jj++){
            wf[ii][jj] = new double[nk];
        }
    }

    for(int ii=0; ii<ni;ii++)
        i[ii] = -imax + di*ii;
    for(int jj=0; jj<nj;jj++){
        j[jj] = -jmax + dj*jj;
    }
    for(int kk=0; kk<nk;kk++)
        k[kk] = -kmax + dk*kk;

    int num_process;
    MPI_Comm_rank(MPI_COMM_WORLD,&num_process);

    for(int ii=0; ii<ni/size_process;ii++){
        for(int jj=0;jj<nj;jj++){
            for(int kk=0;kk<nk;kk++){
                int idx = ii + num_process*ni/size_process;
                wf[ii][jj][kk] = gauss(i[idx],j[jj],k[kk]);
            }
        }
    }

    double integral = 0;
    for(int ii=0; ii<ni/size_process;ii++){
        for(int jj=0;jj<nj;jj++){
            for(int kk=0;kk<nk;kk++){
                integral += wf[ii][jj][kk]*wf[ii][jj][kk]*di*dj*dk;
            }
        }
    }
    MPI_Send(&integral,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
    if(num_process ==0){
        double *integral_arr;
        integral_arr = new double[size_process];
        integral_arr[0] = integral;
        for(int p=1;p<size_process;p++){
           MPI_Recv(&integral_arr[p],1,MPI_DOUBLE,p,0,MPI_COMM_WORLD,&status); 
        }
        double sum = 0;
        for(int p=0;p<size_process;p++){
            sum += integral_arr[p];
        }
        std::cout<<"Integral: "<<sum<<std::endl;
        norm = sum;
        std::cout<<norm<<std::endl;
    }
    if(num_process==0){ 
        for(int p=1;p<size_process;p++){
            MPI_Send(&norm,1,MPI_DOUBLE,p,1,MPI_COMM_WORLD);
        }
    }
    if(num_process!=0)
        MPI_Recv(&norm,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD,&status);

    std::cout<<"n: "<<num_process<<" norm: "<<norm<<std::endl;
    for(int ii=0; ii<ni/size_process;ii++){
        for(int jj=0;jj<nj;jj++){
            for(int kk=0;kk<nk;kk++){
                wf[ii][jj][kk] /= sqrt(norm);
            }
        }
    }
    integral = 0.0;
    for(int ii=0; ii<ni/size_process;ii++){
        for(int jj=0;jj<nj;jj++){
            for(int kk=0;kk<nk;kk++){
                integral += wf[ii][jj][kk]*wf[ii][jj][kk]*di*dj*dk;
            }
        }
    }
    MPI_Send(&integral,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
    if(num_process ==0){
        double *integral_arr;
        integral_arr = new double[size_process];
        integral_arr[0] = integral;
        for(int p=1;p<size_process;p++){
           MPI_Recv(&integral_arr[p],1,MPI_DOUBLE,p,0,MPI_COMM_WORLD,&status); 
        }
        double sum = 0;
        for(int p=0;p<size_process;p++){
            sum += integral_arr[p];
        }
        std::cout<<"Integral: "<<sum<<std::endl;
    }
    std::cout<<"World size: "<<size_process<<" Num: "<<num_process<<std::endl;
    

    MPI_Finalize();
    return 0;
}
