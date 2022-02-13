#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include "debug.h"
#include "wavefunction.h"

void WF::_geom_X(){
    _apply_mask = &WF::apply_mask_X;
}

void WF::apply_mask_X(cdouble *imask, cdouble *jmask, cdouble *kmask){
    for(int i=0; i<_ni; i++){
        _wf[i][0][0] = _wf[i][0][0]*imask[i];
    }
}


