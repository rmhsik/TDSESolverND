#include <iostream>
#include "probe_factory.h"

Probe* ProbeFactory::create(int geometry, std::string def){
    switch(geometry){
        case X:
            std::cout<<"Probe X\n";
            return new ProbeX(def);
        case XZ:
            std::cout<<"Probe XZ\n";
            return new ProbeXZ(def);
            break;
        case RZ:
            std::cout<<"Probe RZ\n";
            return new ProbeRZ(def);
            break;
    }
}
