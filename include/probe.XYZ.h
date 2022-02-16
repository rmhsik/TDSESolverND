#ifndef PROBE_XYZ_H
#define PROBE_XYZ_H

#include "probe.h"

class ProbeXYZ: public Probe{
    protected:
        void _acc_i(const int idx);
        void _acc_j(const int idx);
        void _acc_k(const int idx);
        void _dip_i(const int idx);
        void _dip_j(const int idx);
        void _dip_k(const int idx);
        void _pop(const int idx);
        void _dens(const int idx);
        void _wf_snap(const int idx);
    public:
        ProbeXYZ();
        ProbeXYZ(std::string def);
        virtual ~ProbeXYZ();
};



#endif
