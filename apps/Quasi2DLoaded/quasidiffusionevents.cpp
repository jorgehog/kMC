#include "quasidiffusionevents.h"

#include "quasidiffusion.h"

using namespace kMC;


MovingWall::MovingWall(const double h0, const double EsMax, const double EsInit, const ivec &heighmap):
    KMCEvent("MovingWall", "h0", true, true),
    m_h0(h0),
    m_h(h0),
    m_EsMax(EsMax),
    m_EsInit(EsInit),
    m_r0(r0FromEs(h0, EsMax, EsInit)),
    m_s0(s0FromEs(h0, EsMax, EsInit)),
    m_heighmap(heighmap),
    m_localPressure(heighmap.size()),
    m_expNegOneOverR0(exp(-1.0/m_r0)),
    m_expOneOverR0(exp(1.0/m_r0))
{
    for (uint site = 0; site < m_heighmap.size(); ++site)
    {
        m_localPressure(site) = localPressureEvaluate(site);
    }
}

MovingWall::~MovingWall()
{
    for (uint i = 0; i < m_heighmap.size(); ++i)
    {
        m_pressureAffectedReactions.at(i).clear();
    }

    m_pressureAffectedReactions.clear();
}

void MovingWall::initialize()
{

    m_mPrev = 0;
    m_pressureAffectedReactions.resize(m_heighmap.size());

    uint i;
    for (SoluteParticle *particle : solver()->particles())
    {
        i = particle->x();

        m_mPrev += exp(m_heighmap(i)/m_r0);

        for (Reaction *reaction : particle->reactions())
        {
            if (reaction->isType("QuasiDiffusionReaction"))
            {
                QuasiDiffusionReaction *qdr = static_cast<QuasiDiffusionReaction*>(reaction);

                if (qdr->pressureAffected())
                {
                    m_pressureAffectedReactions[i].push_back(qdr);
                }
            }
        }
    }

    m_mPrev /= m_heighmap.size();

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

    _updatePressureRates();


    for (SoluteParticle* particle : m_affectedParticles)
    {
        for (Reaction *r : particle->reactions())
        {
            if (r->isAllowed())
            {
                r->registerUpdateFlag(QuasiDiffusionReaction::UpdateFlags::CALCULATE);
            }
        }
    }


    if ((cycle() + 1)%10000 == 0)
    {
        remakeUpdatedValues();
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

    m_dh = m_r0*std::log(m/m_mPrev);

    m_h += m_dh;
    m_mPrev = m;

    setValue(m_h);

}

void MovingWall::_updatePressureRates()
{

    double rateChange;

    double expFac = expSmallArg(-m_dh/m_r0);

    for (uint i = 0; i < m_heighmap.size(); ++i)
    {
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
            BADAssSimpleDump(cycle(), i, localPressure(i), expFac, m_heighmap);
        });
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

