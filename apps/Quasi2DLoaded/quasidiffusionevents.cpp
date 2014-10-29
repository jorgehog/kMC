#include "quasidiffusionevents.h"

#include "quasidiffusion.h"

using namespace kMC;


void EqConc::initialize()
{
    m_expFac = exp(Reaction::beta());

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
    setValue(eqConc());
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

    BADAss(localDissolutionRate, <=, 100, "Pressure too high?");

    if (nDissolutionReactions != 0)
    {
        localDissolutionRate /= nDissolutionReactions;
    }

    m_dissolutionRate += dt()*localDissolutionRate;
    m_neighbours += dt()*dependency<NNeighbors>("nNeighbors")->localValue();
    m_totalTime += dt();

    const double &avgNeighbors = m_neighbours/m_totalTime;
    double avgDissolutionRate = m_dissolutionRate/m_totalTime;

    m_eqConc += dt()*avgDissolutionRate/(4 - avgNeighbors);

}


void ConcEquilibriator::initialize()
{
    m_cPrev = solver()->targetConcentration();

    m_counter = 0;

    m_converged = false;

    for (SoluteParticle *particle : solver()->particles())
    {
        for (Reaction *reaction : particle->reactions())
        {
            if (reaction->hasName("Deposition"))
            {
                m_depositionReactions.push_back(static_cast<Deposition*>(reaction));
            }
        }
    }
}

void ConcEquilibriator::execute()
{
    double eqConc = m_eqConcEvent.value();

    if (m_counter >= m_N)
    {
        for (uint i = 1; i < m_N; ++i)
        {
            m_eqConcValues[i - 1] = m_eqConcValues[i];
        }

        m_eqConcValues[m_N - 1] = eqConc;

        double g = flatness();

        if (g < m_gCrit)
        {
            initiateNextConcentrationLevel();
        }


    }
    else
    {
        m_eqConcValues[m_counter] = eqConc;
        m_counter++;
    }

    if (!m_allConcentrations.empty())
    {
        setValue(meanConcentration());
    }

}


void ConcEquilibriator::initiateNextConcentrationLevel()
{
    double cNew = m_eqConcValues[m_N - 1];

    m_allConcentrations.push_back(cNew);
    m_meanConcentration += cNew;

    uint N = 10;

    if (m_allConcentrations.size() > N)
    {

        double cMean = 0;
        for (uint i = 0; i < N; ++i)
        {
            cMean += m_allConcentrations[m_allConcentrations.size() - (i + 1)];
        }

        cMean /= N;

        if (fabs(cMean - m_cPrev) < m_treshold)
        {
            for (double c : m_allConcentrations)
            {
                cout << "c " << c << endl;
            }

            cout << cMean << " " << cNew << " " << m_cPrev << " " << meanConcentration() << endl;

            cNew = cMean;

            terminateLoop("Concentration converged");
        }

    }

    solver()->setTargetConcentration(cNew);
    m_eqConcEvent.restart();

    for (Deposition *depositionReaction : m_depositionReactions)
    {
        depositionReaction->registerUpdateFlag(QuasiDiffusionReaction::UpdateFlags::CALCULATE);
    }

    m_cPrev = cNew;

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
        localValue += exp(Reaction::beta()*m_system.Eb()*m_system.nNeighbors(site));
    }

    m_cumulant = localValue/m_system.size();

    setValue(cumulant());


}

double Cumulant::cumulant() const
{
    return std::log(m_cumulant)/(Reaction::beta()*m_system.Eb());
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

    setValue(m_sum/T());
}
