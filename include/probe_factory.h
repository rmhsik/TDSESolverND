#ifndef PROBE_FACTORY_H
#define PROBE_FACTORY_H

#include "probe.h"
class ProbeFactory{
    public:
        Probe* create(int geometry, std::string def);
};

#endif
