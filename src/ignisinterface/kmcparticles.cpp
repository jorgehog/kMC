#include "kmcparticles.h"

#include "../kmcsolver.h"

#include "../soluteparticle.h"


using namespace kMC;


uint KMCParticles::count() const
{
    return solver()->particles().size();
}

uint KMCParticles::operator ()(const uint n, const uint d) const
{
    return solver()->particle(n)->r(d);
}
