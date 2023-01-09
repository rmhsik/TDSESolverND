#include <iostream>
#include "probe_factory.h"
#include "probe.XZ.h"

Probe* ProbeFactory::create(int geometry, std::string def){
    switch(geometry){
        case XZ:
            std::cout<<"Probe XZ\n";
            return new ProbeXZ(def);
            break;
    }
}
