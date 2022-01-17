#include <iostream>
#include <omp.h>
#include "tdsesolver.h"
#include "parameters.h"

int main(){
    omp_set_num_threads(n_threads);
    Parameters param;
    TDSESolver tdse(param);
    tdse.ipropagate();
    tdse.propagate();
}
