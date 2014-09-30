#include "quasidiffusionevents.h"

#include "quasidiffusion.h"

using namespace kMC;


void EquilibriumConcentrationEstimator::initialize()
{
    m_expFac = exp(Reaction::beta());

    restart();

    for (SoluteParticle *particle : solver()->particles())
    {
        for (const Reaction *reaction : particle->reactions())
        {
            if (reaction->hasName("Dissolution"))
            {
                m_dissolutionReactions.push_back(static_cast<const Dissolution*>(reaction));
            }
        }
    }
}

void EquilibriumConcentrationEstimator::execute()
{
    if (cycle() != 0)
    {
        setValue(m_expFac*m_eqConc/cycle());
    }
}

void EquilibriumConcentrationEstimator::reset()
{
    double localNeighbors = 0;
    double localDissolutionRate = 0;

    uint nDissolutionReactions = 0;
    for (const Dissolution *dissolutionReaction : m_dissolutionReactions)
    {
        if (dissolutionReaction->isAllowed())
        {
            localDissolutionRate += dissolutionReaction->rate();
            nDissolutionReactions++;
        }

        localNeighbors += dissolutionReaction->nNeighbors();
    }

    if (localDissolutionRate > 100)
    {
        solver()->exit("Pressure too high?");
    }

    localNeighbors /= m_dissolutionReactions.size();
    if (nDissolutionReactions != 0)
    {
        localDissolutionRate /= nDissolutionReactions;
    }

    m_neighbours += localNeighbors;
    m_dissolutionRate += localDissolutionRate;

    double avgNeighbours = m_neighbours/(cycle() + 1);
    double avgDissolutionRate = m_dissolutionRate/(cycle() + 1);

    m_eqConc += avgDissolutionRate/(4 - avgNeighbours);
}
