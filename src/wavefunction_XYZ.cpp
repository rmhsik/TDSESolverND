#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include "debug.h"
#include "wavefunction.h"

void WF::_geom_XYZ(){
    _apply_mask = &WF::apply_mask_XYZ;
}

void WF::apply_mask_XYZ(cdouble *imask, cdouble *jmask, cdouble *kmask){
    for(int i=0;i<_ni;i++){
        for(int j=0;j<_nj;j++){
            for(int k=0;k<_nk;k++){
                _wf[i][j][k] *= imask[i]*jmask[j]*kmask[k];
            }
        }
    }
}

void WF::apply_mask_buf_XYZ(cdouble *imask, cdouble *jmask, cdouble *kmask, const int idx){
    for(int i=0; i<_ni; i++){
        for(int j=0; j<_nj; j++){
            for(int k=0; k<_nk; k++){
                _wf_buf[idx][i][j][k] *= imask[i]*jmask[j]*kmask[k];
            }
        }
    }

}
