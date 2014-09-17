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
    m_localPressure(heighmap.size())
{
    for (uint site = 0; site < m_heighmap.size(); ++site)
    {
        m_localPressure(site) = localPressureEvaluate(site);
    }
}

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

    for (uint i = 0; i < m_heighmap.size(); ++i)
    {
        m_localPressure(i) = localPressureEvaluate(i);
    }

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

    for (QuasiDiffusionReaction* r : m_pressureAffectedReactions)
    {
        BADAssBool(r->pressureAffected());

        if (!r->isAllowed())
        {
            continue;
        }

        else if (isAffected(r->reactant()))
        {
            continue;
        }

        double R0 = r->rate()/exp(-Reaction::beta()*localPressureOld(r->x()));

        r->changeRate([&] (const double rate)
        {
            double delta = m_localPressure(r->x()) - localPressureOld(r->x());
            double arg = -Reaction::beta()*delta;
            double arg2 = arg*arg;
            double approx = 1 + arg + 0.5*arg2 + 1./6*arg*arg2;

            BADAssClose(exp(arg), approx, 1E-3);

            return rate*approx;
        });

        double R2 = r->rate();
        r->registerUpdateFlag(QuasiDiffusionReaction::UpdateFlags::CALCULATE);
        double R = r->calcRate();

        r->forceNewRate(R);

        //DETTE SKJER KUN I SYKEL NULL. HVA ER GALT? NOE SOM IKKE ER SOM DET SKAL.
        if (fabs(R - R2) > 0.001)
        {
            cout << R << " " << R2 << endl;
            cout << r->info() << endl;
            cout << solver()->solverEvent()->selectedReaction()->info() << endl;
            cout << solver()->solverEvent()->cycle() << endl;
            cout << "----" << endl;
        }

        double R3 = R/exp(-Reaction::beta()*localPressure(r->x()));

        double n0 = (-log(R0) - 1)/r->Eb();
        double n3 = (-log(R3) - 1)/r->Eb();

        BADAssClose(n0, n3, 1E-10, "fail", [&] ()
        {
            cout << solver()->solverEvent()->selectedReaction()->info() << endl;
            cout << "---" << endl;
            cout << r->info() << endl;

            BADAssSimpleDump(solver()->solverEvent()->cycle(),
                             r->nNeighbors(),
                             n0,
                             n3,
                             localPressure(r->x()),
                             localPressureOld(r->x()));
        });
        BADAssClose(R, R2, 1E-3);

    }
}

