#pragma once

#include "kmcparticles.h"

#include <ignis/include/ignis.h>

using ignis::Event;

namespace kMC
{

class DiffusionReaction;

class KMCEvent : public Event<uint>
{
    REGISTER_POSITIONHANDLER(KMCParticles, uint)

public:

    KMCEvent(std::string type = "kMCEvent",
             std::string unit = "",
             bool doOutput = false,
             bool toFile = false) :
        Event<uint>(type,
                    unit,
                    doOutput,
                    toFile)
    {

    }

    virtual ~KMCEvent() {}

protected:

    KMCSolver *solver() const
    {
        return registeredHandler().solver();
    }

    const DiffusionReaction *lastReaction() const;

    const uint &NX() const;

    const uint &NY() const;

    const uint &NZ() const;

};

}
