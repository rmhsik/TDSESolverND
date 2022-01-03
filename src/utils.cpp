#include "utils.h"
#include "debug.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <stdlib.h>

int calc_n_elem(const std::string &path){
    std::ifstream file;
    file.open(path);
    std::string line;
    int n_elem = 0;
    if(file.is_open()){
        while(getline(file,line)){
            n_elem++;
        }
        file.close();
   } 
    else{debug0("[calc_n_elem] Error opening file.\n"); exit(1);} 
   return n_elem;
} 


std::string define_filepath(std::string &base){
    int i = 0; 
    //std::string basepath = "results/prop_acc";
    std::string ext = ".dat";
    std::string path = base + std::to_string(i) + ext;
    std::ifstream file; 
    file.open(path);
    while(file.is_open()){
        i++;
        file.close();
        path = base + std::to_string(i) +ext;
        file.open(path);
    }
    file.close();
    return path;
}

template <class T>
void write_array(T *arr, const int n, const std::string &path){
    std::ofstream outfile;
    outfile.open(path);
    if(outfile.is_open()){    
        for(int i=0; i<n; i++){
            std::ostringstream doubleStr;
            doubleStr<<std::fixed<<std::setprecision(12);
            doubleStr<<arr[i];
            outfile<<doubleStr.str()<<std::endl;
        }
    }
    else{debug0("[write_vector] Error opening file.\n"); exit(1);}
}


template <class T>
std::tuple<T* ,double> linspace(T xi, T xf, int n){
    double dx = abs(xf -xi)/(double)(n-1); 
    T* grid = new T[n];
    if (n==0){}
    else if(n==1){ grid[0]==xi;}
    else{
        for(int i=0; i<n-1;i++){
            grid[i] = xi + i*dx;
        }
        grid[n-1] = xf;
    }
    return std::make_tuple(grid, dx);
}

template void write_array<int>(int *arr,const int n, const std::string &path);
template void write_array<double>(double *arr, const int n, const std::string &path);
template void write_array<std::complex<double>>(std::complex<double> *arr,const int n, const std::string &path);

template std::tuple<int*,double> linspace<int>(int, int, int);
template std::tuple<double*,double> linspace<double>(double, double, int);
template std::tuple<std::complex<double>*,double> linspace<std::complex<double>>(std::complex<double>, std::complex<double>, int);
