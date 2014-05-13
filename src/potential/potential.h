#pragma once

#include <sys/types.h>

namespace kMC
{

class SoluteParticle;

class Potential
{
public:

    Potential()
    {

    }

    virtual ~Potential() {}

    virtual void initialize() {}

    virtual double valueAt(const double x, const double y, const double z) = 0;

    virtual double evaluateFor(SoluteParticle *particle) = 0;

    virtual double evaluateSaddleFor(SoluteParticle *particle, const uint dx, const uint dy, const uint dz) = 0;

    virtual double onNeighborChange(SoluteParticle *neighbor, const uint dx, const uint dy, const uint dz) = 0;

};

}
