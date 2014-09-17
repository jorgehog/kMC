#include "quasidiffusionevents.h"

#include "quasidiffusion.h"

using namespace kMC;


void MovingWall::initialize()
{
    for (SoluteParticle *particle : solver()->particles())
    {
        for (Reaction *reaction : particle->reactions())
        {
            if (reaction->isType("QuasiDiffusionReaction"))
            {
                QuasiDiffusionReaction *qdr = static_cast<QuasiDiffusionReaction*>(reaction);

                if (qdr->pressureAffected())
                {
                    m_pressureAffectedReactions.push_back(qdr);
                }
            }
        }
    }

    for (uint site = 0; site < m_heighmap.size(); ++site)
    {
        m_localPressure(site) = localPressureEvaluate(site);
    }
}

void MovingWall::execute()
{
    for (uint i = 0; i < m_heighmap.size(); ++i)
    {
        if (!isAffected(solver()->particle(i)))
        {
            BADAssClose(localPressureEvaluate(i), localPressure(i), 1E-3);
        }
    }

    _rescaleHeight();

    _updatePressureRates();

}

void MovingWall::reset()
{
    for (SoluteParticle *particle : SoluteParticle::affectedParticles())
    {
        BADAssClose(localPressureEvaluate(particle->x()), localPressure(particle->x()), 1E-3);
    }

    for (QuasiDiffusionReaction *r : m_pressureAffectedReactions)
    {
        if (r->reactant()->isAffected() || !r->isAllowed())
        {
            continue;
        }

        BADAssClose(r->calcRateBruteForce(), r->rate(), 1E-3);
        BADAssClose(r->calcRateBruteForce(), r->calcRate(), 1E-3);
    }

}

void MovingWall::_updatePressureRates()
{
    for (SoluteParticle *particle : solver()->particles())
    {
        m_localPressure(particle->x()) = localPressureEvaluate(particle->x());
    }

    for (QuasiDiffusionReaction* r : m_pressureAffectedReactions)
    {
        if (!r->isAllowed() || isAffected(r->reactant()))
        {
            continue;
        }

        r->setRate();
        r->registerUpdateFlag(QuasiDiffusionReaction::UpdateFlags::SKIP);
    }

    m_affectedParticles.clear();
}

