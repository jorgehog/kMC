#include "quasidiffusionsystem.h"

#include "quasidiffusion.h"

using namespace kMC;


void QuasiDiffusionSystem::setHeights(const ivec &heights)
{
    m_heights = heights;

    m_wallEvent.onHeightsSet();

    for (Reaction *reaction : KMCSolver::instance()->allPossibleReactions())
    {
        reaction->registerUpdateFlag(QuasiDiffusionReaction::UpdateFlags::CALCULATE);
    }

}
