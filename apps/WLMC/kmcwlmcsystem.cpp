#include "kmcwlmcsystem.h"

#include <kMC>

using namespace kMC;


KMCWLMCSystem::KMCWLMCSystem(KMCSolver *solver,
                             const uint movesPerSampling,
                             const double flatnessCriterion,
                             const uint overlap,
                             const uint minWindowSize,
                             const uint windowIncrementSize,
                             const double *f) :
    WLMCSystem(SoluteParticle::nParticles(),
               solver->NX(),
               solver->NY(),
               solver->NZ(),
               movesPerSampling,
               flatnessCriterion,
               overlap,
               minWindowSize,
               windowIncrementSize,
               f,
               [] () {return KMC_RNG_UNIFORM();}),
    m_solver(solver)
{

}

bool KMCWLMCSystem::isOccupiedLoction(const uint x, const uint y, const uint z) const
{
    return m_solver->getSite(x, y, z)->isActive();
}

double KMCWLMCSystem::getValue(const uint particleIndex) const
{
    return m_solver->particle(particleIndex)->energy();
}

double KMCWLMCSystem::getTotalValue() const
{
    return SoluteParticle::totalEnergy();
}

double KMCWLMCSystem::getValueDifference(const uint particleIndex, const uint xd, const uint yd, const uint zd) const
{
    double eNew = 0;

    SoluteParticle *particle = m_solver->particle(particleIndex);

    Site::forEachNeighborDo_sendIndices(xd, yd, zd, [&particle, &eNew] (Site *neighbor, uint i, uint j, uint k)
    {
        if (neighbor->isActive())
        {
            if (neighbor->associatedParticle() != particle)
            {
                eNew += DiffusionReaction::potential(i, j, k);
            }
        }
    });

    return 2*(eNew - particle->energy());
}

void KMCWLMCSystem::changePosition(const uint particleIndex, const uint xd, const uint yd, const uint zd)
{
    m_solver->particle(particleIndex)->changePosition(xd, yd, zd);
}
