#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include "debug.h"
#include "wavefunction.h"

void WF::_geom_XZ(){
    _apply_mask = &WF::apply_mask_XZ;
}

void WF::apply_mask_XZ(cdouble *imask, cdouble *kmask){
    for(int i=0; i<_ni; i++){
        for(int k=0; k<_nk;k++)
            _wf[i][k] = _wf[i][k]*imask[i]*kmask[k];
    }
}

void WF::apply_mask_buf_XZ(cdouble *imask, cdouble *kmask,const int idx){
    for(int i=0; i<_ni; i++){
        for(int k=0; k<_nk;k++)
            _wf_buf[idx][i][k] *= imask[i]*kmask[k];
    }
}


