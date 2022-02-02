#ifndef PROBE_RZ_H
#define PROBE_RZ_H

#include "probe.h"

class ProbeRZ: public Probe{
    protected:
        void _acc_i(const int idx);
        void _acc_k(const int idx);
        void _dip_i(const int idx);
        void _dip_k(const int idx);
        void _pop(const int idx);
        void _dens(const int idx);
    public:
        ProbeRZ();
        ProbeRZ(std::string def);
};

#endif
