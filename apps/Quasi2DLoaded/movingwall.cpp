#include "quasidiffusionevents.h"

#include "quasidiffusion.h"

using namespace kMC;


MovingWall::MovingWall(const double E0,
                       const double sigma0,
                       const double r0,
                       const ivec &heighmap):
    KMCEvent("MovingWall", "h0", true, true),
    m_r0(r0),
    m_s0(sigma0),
    m_E0(E0),
    m_heighmap(heighmap),
    m_localPressure(heighmap.size(), fill::zeros)
{

}

MovingWall::~MovingWall()
{
    for (uint i = 0; i < m_pressureAffectedReactions.size(); ++i)
    {
        m_pressureAffectedReactions.at(i).clear();
    }

    m_pressureAffectedReactions.clear();

}

void MovingWall::initialize()
{
    BADAssBool(isActive() && hasStarted());

    m_mPrev = 0;
    m_pressureAffectedReactions.resize(m_heighmap.size());

    uint site;
    for (SoluteParticle *particle : solver()->particles())
    {
        site = particle->x();

        m_mPrev += exp(m_heighmap(site)/m_r0);
        m_pressureAffectedReactions[site].clear();

        for (Reaction *reaction : particle->reactions())
        {
            if (reaction->isType("QuasiDiffusionReaction"))
            {
                QuasiDiffusionReaction *qdr = static_cast<QuasiDiffusionReaction*>(reaction);

                if (qdr->pressureAffected())
                {
                    m_pressureAffectedReactions[site].push_back(qdr);
                }
            }
        }
    }

    m_mPrev /= m_heighmap.size();

    m_h = m_r0*std::log(m_heighmap.size()*m_s0*m_mPrev/(-m_E0));
    cout << m_h << endl;

    for (uint site = 0; site < m_heighmap.size(); ++site)
    {
        m_localPressure(site) = localPressureEvaluate(site);

        for (Reaction *reaction : solver()->particle(site)->reactions())
        {
            reaction->registerUpdateFlag(QuasiDiffusionReaction::UpdateFlags::CALCULATE);
        }
    }

    solver()->getRateVariables();

    BADAssClose(pressureEnergySum(), m_E0, 1E-5);

}

void MovingWall::reset()
{

    //Check that everything is updated correctly from previous runs.
#ifndef NDEBUG
    for (uint i = 0; i < m_heighmap.size(); ++i)
    {
        if (!isAffected(solver()->particle(i)))
        {
            BADAssClose(localPressureEvaluate(i), localPressure(i), 1E-3);
        }
    }
#endif

    //Calculate new height of wall to conserve total force.
    _rescaleHeight();

    //Check if the move allowed any new reactions to occur.
    _locateNewAffected();

    //Check if the move blocked any active reactions.
    _removeBlocked();

    _updatePressureRates();

    BADAssClose(pressureEnergySum(), m_E0, 1E-5);


    if ((cycle() + 1)%10000 == 0)
    {
        remakeUpdatedValues();
    }


    for (SoluteParticle *particle : m_affectedParticles)
    {
        for (Reaction *reaction : particle->reactions())
        {
            if (reaction->isAllowed())
            {
                reaction->registerUpdateFlag(QuasiDiffusionReaction::UpdateFlags::CALCULATE);
                reaction->setRate();
            }
        }
    }

    //if any of the reactions are enabled and disabled, we need to remake the event list.
    //This should really be automated in the KMCSolver
    if (m_changed != 0)
    {
        KMCSolver::instance()->reshuffleReactions();
        KMCSolver::instance()->remakeAccuAllRates();
    }

    m_affectedParticles.clear();
}

void MovingWall::_rescaleHeight()
{

    double m = 0;
    for (uint i = 0; i < m_heighmap.size(); ++i)
    {
        m += exp(m_heighmap(i)/m_r0);
    }

    m /= m_heighmap.size();

    m_dh = m_r0*std::log(m/m_mPrev);

    m_mPrev = m;

    m_h += m_dh;

    setValue(m_h - dependency("height")->value());

}

void MovingWall::_updatePressureRates()
{

    double rateChange;

    double expFac = expSmallArg(-m_dh/m_r0);

    for (uint i = 0; i < m_heighmap.size(); ++i)
    {

        if (!isAffected(solver()->particle(i)))
        {

            rateChange = expSmallArg(-m_system->alpha()*m_localPressure(i)*(expFac - 1));

            //For every affected particle we update only those who include the pressure term.
            //Vector is set up in initialize based on virtual reaction function isPressureAffected().
            for (auto &r : m_pressureAffectedReactions[i])
            {
                BADAssBool(r->pressureAffected(), "This reaction does not belong here.", [&r] ()
                {
                    cout << r->info() << endl;
                });

                if (!r->isAllowed())
                {
                    continue;
                }

                r->changeRate(r->rate()*rateChange);

            }
            m_localPressure(i) *= expFac;
        }
        else
        {
            m_localPressure(i) = localPressureEvaluate(i);
        }

        BADAssClose(localPressureEvaluate(i), m_localPressure(i), 1E-3, "incorrect pressure update", [&] ()
        {
            BADAssSimpleDump(cycle(), i, localPressure(i), expFac, m_dh);
        });
    }

}

void MovingWall::_locateNewAffected()
{
    m_changed = 0;

    for (SoluteParticle *particle : solver()->particles())
    {
        for (Reaction *reaction : particle->reactions())
        {
            if (reaction->isAllowed() && (reaction->rate() == Reaction::UNSET_RATE))
            {
                m_changed++;
                m_affectedParticles.insert(particle);
                break;
            }
        }
    }
}

void MovingWall::_removeBlocked()
{
    for (SoluteParticle *particle : solver()->particles())
    {
        for (Reaction *reaction : particle->reactions())
        {
            if (!reaction->isAllowed() && (reaction->rate() != Reaction::UNSET_RATE))
            {
                m_changed++;
                reaction->disable();
                reaction->resetUpdateFlag();
            }
        }
    }

}

void MovingWall::remakeUpdatedValues()
{

    recalculateAllPressures();

    for (Reaction *reaction : solver()->allPossibleReactions())
    {
        reaction->registerUpdateFlag(QuasiDiffusionReaction::UpdateFlags::CALCULATE);
        reaction->setRate();
    }
}

void MovingWall::recalculateAllPressures()
{
    for (uint i = 0; i < m_heighmap.size(); ++i)
    {
        m_localPressure(i) = localPressureEvaluate(i);
    }
}

