#include "probe_factory.h"

Probe* ProbeFactory::create(int geometry, std::string def){
    switch(geometry){
        case XZ:
            return new ProbeXZ(def);
            break;
        case RZ:
            return new ProbeRZ(def);
            break;
    }
}
