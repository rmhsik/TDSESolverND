#ifndef PROBE_X_H
#define PROBE_X_H

#include "probe.h"

class ProbeX: public Probe{
    protected:
        void _acc_i(const int idx);
        void _acc_k(const int idx);
        void _dip_i(const int idx);
        void _dip_k(const int idx);
        void _pop(const int idx);
        void _dens(const int idx);
    public:
        ProbeX();
        ProbeX(std::string def);
};

#endif
