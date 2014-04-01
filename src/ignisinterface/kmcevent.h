#pragma once

#include "kmcparticles.h"

#include "../../ignis/include/ignis.h"

using ignis::Event;

namespace kMC
{

class KMCParticles;

class KMCEvent : public Event<uint>
{
    REGISTER_POSITIONHANDLER(KMCParticles, uint)

public:

    KMCEvent(string type = "kMCEvent",
               string unit = "",
               bool doOutput = false,
               bool toFile = false) :
        Event<uint>(type,
                    unit,
                    doOutput,
                    toFile)
    {

    }

protected:

    static KMCSolver *solver()
    {
        return particles().solver();
    }

};

}
