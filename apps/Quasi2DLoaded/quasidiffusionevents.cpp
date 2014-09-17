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
    for (SoluteParticle *particle : solver()->particles())
    {
        BADAssClose(localPressureEvaluate(particle->x()), localPressure(particle->x()), 1E-3);
    }

    m_affectedParticles.clear();
}

void MovingWall::_updatePressureRates()
{
    vec localPressureOld = m_localPressure;

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

        r->changeRate([&] (const double rate)
        {
            return rate*exp(-Reaction::beta()*(m_localPressure(r->reactant()->x()) - localPressureOld(r->reactant()->x())));
        });

        r->registerUpdateFlag(QuasiDiffusionReaction::UpdateFlags::CALCULATE);
    }
}

