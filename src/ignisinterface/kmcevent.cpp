#include "kmcevent.h"

#include "../kmcsolver.h"

#include "solverevent.h"

using namespace kMC;


const double &KMCEvent::dt() const
{
    return solver()->solverEvent()->lastTimeStep();
}

const double &KMCEvent::T() const
{
    return solver()->solverEvent()->totalTime();
}

const DiffusionReaction *KMCEvent::lastReaction() const
{
    return reinterpret_cast<const DiffusionReaction*>(solver()->solverEvent()->selectedReaction());
}

const uint &KMCEvent::NX() const
{
    return solver()->NX();
}

const uint &KMCEvent::NY() const
{
    return solver()->NY();
}

const uint &KMCEvent::NZ() const
{
    return solver()->NZ();
}
