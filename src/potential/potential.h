#pragma once

#include <sys/types.h>

namespace kMC
{

class SoluteParticle;
class DiffusionReaction;

class Potential
{
public:

    Potential()
    {

    }

    virtual ~Potential() {}

    virtual void initialize() {}

    virtual double valueAt(const double x, const double y, const double z) = 0;

    virtual double evaluateFor(const SoluteParticle *particle) = 0;

    virtual double evaluateSaddleFor(const DiffusionReaction *currentReaction) = 0;

    virtual double onNeighborChange(SoluteParticle *particle,
                                    const SoluteParticle *neighbor,
                                    const uint dx,
                                    const uint dy,
                                    const uint dz,
                                    int sign) = 0;

};

}
