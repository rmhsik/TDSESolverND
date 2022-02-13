#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include "debug.h"
#include "wavefunction.h"

void WF::_geom_RZ(){
    _apply_mask = &WF::apply_mask_RZ;
}

void WF::apply_mask_RZ(cdouble *imask, cdouble *jmask, cdouble *kmask){
    for(int i=0; i<_ni; i++){
        for(int k=0; k<_nk;k++)
            _wf[i][0][k] = _wf[i][0][k]*imask[i]*kmask[k];
    }
}


void WF::apply_mask_buf_RZ(cdouble *imask, cdouble *jmask, cdouble *kmask,const int idx){
    for(int i=0; i<_ni; i++){
        for(int k=0; k<_nk;k++)
            _wf_buf[idx][i][0][k] *= imask[i]*kmask[k];
    }
}

