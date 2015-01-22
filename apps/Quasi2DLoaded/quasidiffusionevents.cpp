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

    m_totalTime += dt();

    m_dissolutionRate += dt()*localDissolutionRate;

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

    m_doAverage = false;

    m_averageMu = 0;

    m_averageMu2Sum = 0;

    m_averageMuCount = 0;

    m_finalized = false;

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

    m_counter++;
    if (m_counter >= m_N)
    {
        initiateNextConcentrationLevel(shift);
        m_counter = 0;
    }

    setValue(m_system.mu());

}


void ConcEquilibriator::initiateNextConcentrationLevel(const double shift)
{

    double newMu = m_system.mu() + shift;

    //Cannot use 'else' because this could occur even when the prev test goes through
    if (m_doAverage)
    {
        m_averageMu += newMu;
        m_averageMu2Sum += newMu*newMu;
        m_averageMuCount++;
    }

    m_shifts.push_back(shift);
    m_values.push_back(m_system.mu());

    m_system.setMu(newMu);

    conv_to<vec>::from(m_shifts).eval().save("/tmp/shifts.arma");
    conv_to<vec>::from(m_values).eval().save("/tmp/values.arma");


    if (m_averageMuCount == m_nRounds)
    {
        terminateLoop("Concentration converged");

        finalizeAverages();
    }

    if (!m_doAverage)
    {
        if (m_prevShift != 0 && shift != 0)
        {
            if (m_prevShift/shift < 0)
            {
                m_doAverage = true;
            }
        }
    }

    m_prevShift = shift;

    m_eqConcEvent.restart();

    for (QuasiDiffusionReaction *reaction : m_affectedReactions)
    {
        reaction->registerUpdateFlag(QuasiDiffusionReaction::UpdateFlags::CALCULATE);
    }
}

void ConcEquilibriator::finalizeAverages()
{
    if (m_finalized)
    {
        return;
    }

    for (uint i = 0; i < m_shifts.size(); ++i)
    {
        cout << m_shifts.at(i) << " " << m_values.at(i) << endl;
    }

    const uint &N = m_averageMuCount;

    cout << "N = " << N << endl;

    m_averageMu /= N;

    m_error = sqrt(1.0/(N - 1)*(m_averageMu2Sum - m_averageMu*m_averageMu*N));

    cout << "Average = " << m_averageMu << " error = " << m_error << endl;

    m_finalized = true;
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
