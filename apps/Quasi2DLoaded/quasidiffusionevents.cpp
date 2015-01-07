#include "quasidiffusionevents.h"

#include "quasidiffusion.h"

using namespace kMC;


void EqConc::initialize()
{
    resetCounters();

    m_dissolutionReactions.clear();
    for (SoluteParticle *particle : solver()->particles())
    {
        for (Reaction *reaction : particle->reactions())
        {
            if (reaction->hasName("Dissolution"))
            {
                m_dissolutionReactions.push_back(static_cast<Dissolution*>(reaction));
            }
        }
    }

    update();
}

void EqConc::execute()
{
    setValue(dMu());
}

void EqConc::reset()
{
    update();
}

void EqConc::restart()
{
    resetCounters();
}

void EqConc::update()
{
    double localDissolutionRate = 0;

    uint nDissolutionReactions = 0;
    for (const Dissolution *dissolutionReaction : m_dissolutionReactions)
    {
        if (dissolutionReaction->isAllowed())
        {
            localDissolutionRate += dissolutionReaction->rate();
            nDissolutionReactions++;
        }
    }

    if (nDissolutionReactions != 0)
    {
        localDissolutionRate /= nDissolutionReactions;
    }

    m_dissolutionRate += dt()*localDissolutionRate;
    m_totalTime += dt();

    double avgDissolutionRate = m_dissolutionRate/m_totalTime;

    double inverseKStar;

    if (m_shadowing)
    {
        m_neighbours += dt()*dependency<NNeighbors>("nNeighbors")->localValue();
        const double &avgNeighbors = m_neighbours/m_totalTime;

        inverseKStar = 1./DepositionMirrorImageArhenius::promoteFactor(avgNeighbors);
    }
    else
    {
        inverseKStar = 1./DepositionMirrorImageArheniusNoShadowing::promoteFactor();
    }

    m_dMu = avgDissolutionRate*inverseKStar;

}


void ConcEquilibriator::initialize()
{
    m_counter = 0;

    m_initialHeights = m_system.heights();

    m_prevShift = 0;
    m_nswaps = 0;

    for (SoluteParticle *particle : solver()->particles())
    {
        for (Reaction *reaction : particle->reactions())
        {
            if (!reaction->hasName("Deposition"))
            {
                m_affectedReactions.push_back(static_cast<QuasiDiffusionReaction*>(reaction));
            }
        }
    }
}

void ConcEquilibriator::execute()
{
    double shift = m_eqConcEvent.value();

    if (m_counter >= m_N)
    {
        for (uint i = 1; i < m_N; ++i)
        {
            m_logCShiftValues[i - 1] = m_logCShiftValues[i];
        }

        m_logCShiftValues[m_N - 1] = shift;

        double g = flatness();

        if (g < m_gCrit)
        {
            initiateNextConcentrationLevel();
            m_counter = 0;
        }
    }
    else
    {
        m_logCShiftValues[m_counter] = shift;
        m_counter++;
    }

    setValue(m_system.mu());

}


void ConcEquilibriator::initiateNextConcentrationLevel()
{
    double shift = m_logCShiftValues[m_N - 1];
    uint swapmax = 100;

    if (m_prevShift != 0 && shift != 0)
    {
        if (m_prevShift/shift < 0)
        {
            if (m_nswaps < swapmax)
            {
                m_nswaps++;
            }
        }
        else
        {
            if (m_nswaps != 0)
            {
                m_nswaps--;
            }
        }
    }

    m_prevShift = shift;

    shift /= (m_nswaps + 1);

    m_system.setMu(m_system.mu() + shift);

    m_shifts.push_back(shift);
    m_values.push_back(m_system.mu());

    conv_to<vec>::from(m_shifts).eval().save("/tmp/shifts.arma");
    conv_to<vec>::from(m_values).eval().save("/tmp/values.arma");

    if (fabs(shift) < m_treshold)
    {
        terminateLoop("Concentration converged");

        for (uint i = 0; i < m_shifts.size(); ++i)
        {
            cout << m_shifts.at(i) << " " << m_values.at(i) << endl;
        }
    }


//    m_system.setHeights(m_initialHeights);
    m_eqConcEvent.restart();

    for (QuasiDiffusionReaction *reaction : m_affectedReactions)
    {
        reaction->registerUpdateFlag(QuasiDiffusionReaction::UpdateFlags::CALCULATE);
    }
}


void NNeighbors::execute()
{
    m_localValue = 0;

    for (uint site = 0; site < m_system.size(); ++site)
    {
        m_localValue += m_system.nNeighbors(site);
    }

    m_localValue /= m_system.size();

    m_sum += m_localValue;

    setValue(m_sum/(cycle() + 1));
}


void Cumulant::execute()
{
    double localValue = 0;

    for (uint site = 0; site < m_system.size(); ++site)
    {
        localValue += exp(m_system.alpha()*m_system.nNeighbors(site));
    }

    m_cumulant = localValue/m_system.size();

    setValue(cumulant());
}

double Cumulant::cumulant() const
{
    return std::log(m_cumulant)/m_system.alpha();
}


void SurfaceSize::execute()
{
    m_localValue = 0;

    for (uint site = 0; site < m_system.size(); ++site)
    {
        m_localValue += 0.5*abs(m_system.heights(site) - m_system.heights(m_system.rightSite(site, 1)));
    }

    m_localValue /= m_system.size();

    m_sum += dt()*m_localValue;

    setValue(m_sum/(T() - m_T0));
}
