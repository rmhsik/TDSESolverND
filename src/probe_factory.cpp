#include <iostream>
#include "probe_factory.h"
#include "probe.XZ.h"

Probe* ProbeFactory::create(int geometry, std::string def){
    switch(geometry){
        case XZ:
            return new ProbeXZ(def);
            break;
    }
}
