#include "quasidiffusionevents.h"

#include "quasidiffusion.h"

using namespace kMC;


MovingWall::MovingWall(const double h0, const double EsMax, const double EsInit, const ivec &heighmap):
    KMCEvent("MovingWall", "h0", true, true),
    m_inContact(true),
    m_EsMax(EsMax),
    m_EsInit(EsInit),
    m_r0(r0FromEs(h0, EsMax, EsInit)),
    m_s0(s0FromEs(h0, EsMax, EsInit)),
    m_E0(heighmap.size()*_pressureExpression(h0)),
    m_h0(h0),
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

    if (m_h <= m_heighmap.max())
    {
        m_h = m_heighmap.max();
        m_inContact = true;
    }
    else
    {
        m_inContact = false;
    }

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

void MovingWall::execute()
{
#ifndef NDEBUG
    for (uint i = 0; i < m_heighmap.size(); ++i)
    {
        if (!isAffected(solver()->particle(i)))
        {
            BADAssClose(localPressureEvaluate(i), localPressure(i), 1E-3);
        }
    }
#endif

    _rescaleHeight();

    _locateNewAffected();

    _updatePressureRates();

    if (!m_inContact)
    {
        BADAssClose(pressureEnergySum(), m_E0, 1E-5);
    }

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
            }
        }
    }
}

void MovingWall::reset()
{
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

    double trialHeight;

    if (m_inContact)
    {
        trialHeight = m_r0*std::log(m_heighmap.size()*m_s0*m/(-m_E0));
    }

    else
    {
        m_dh = m_r0*std::log(m/m_mPrev);

        trialHeight = m_h + m_dh;
    }


    if (trialHeight <= m_heighmap.max())
    {

#ifndef NDEBUG
        cout << "CONTACT: Unable to shift wall. " << m_h << " / " << m_heighmap.max() << endl;
#endif

        m_dh = m_heighmap.max() - m_h;
        m_h += m_dh;

        if (!m_inContact)
        {

            m_inContact = true;
            remakeUpdatedValues();

            m_dh = 0;
        }

    }
    else
    {
        m_h = trialHeight;

        if (m_inContact)
        {
            m_inContact = false;
            remakeUpdatedValues();

            m_dh = 0;
        }

        m_mPrev = m;
    }

    setValue(m_h - dependency("height")->value());

}

void MovingWall::_updatePressureRates()
{

    double rateChange;

    double expFac = expSmallArg(-m_dh/m_r0);

    for (uint i = 0; i < m_heighmap.size(); ++i)
    {
        BADAss(m_heighmap(i), <=, m_h);

        if (!isAffected(solver()->particle(i)))
        {
            rateChange = expSmallArg(-Reaction::beta()*m_localPressure(i)*(expFac - 1));

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
    for (SoluteParticle *particle : solver()->particles())
    {
        for (Reaction *reaction : particle->reactions())
        {
            if (reaction->isAllowed() && (reaction->rate() == Reaction::UNSET_RATE))
            {
                m_affectedParticles.insert(particle);
                break;
            }
        }
    }
}

void MovingWall::remakeUpdatedValues()
{

    for (uint i = 0; i < m_heighmap.size(); ++i)
    {
        m_localPressure(i) = localPressureEvaluate(i);
    }

    for (Reaction *reaction : solver()->allPossibleReactions())
    {
        reaction->registerUpdateFlag(QuasiDiffusionReaction::UpdateFlags::CALCULATE);
        reaction->setRate();
    }
}

