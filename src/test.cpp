#include <iostream>
#include <omp.h>
#include "tdsesolver.h"
#include "parameters.h"

int main(){
    Parameters param;
    param.print();
    omp_set_num_threads(param.n_threads);
    TDSESolver tdse(param);
    tdse.ipropagate();
    tdse.propagate();
}
