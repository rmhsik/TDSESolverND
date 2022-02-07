#ifndef PROBE_XZ_H
#define PROBE_XZ_H

#include "probe.h"

class ProbeXZ: public Probe{
    protected:
        void _acc_i(const int idx);
        void _acc_k(const int idx);
        void _dip_i(const int idx);
        void _dip_k(const int idx);
        void _pop(const int idx);
        void _dens(const int idx);
        void _wf_snap(const int idx);
    public:
        ProbeXZ();
        ProbeXZ(std::string def);
};

#endif
