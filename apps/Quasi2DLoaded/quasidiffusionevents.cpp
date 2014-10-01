#include "quasidiffusionevents.h"

#include "quasidiffusion.h"

using namespace kMC;


void EqConc::initialize()
{
    m_expFac = exp(Reaction::beta());

    resetCounters();

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
    double localNeighbors = 0;
    double localDissolutionRate = 0;

    uint nDissolutionReactions = 0;
    for (const Dissolution *dissolutionReaction : m_dissolutionReactions)
    {
        if (dissolutionReaction->isAllowed())
        {
            localDissolutionRate += dissolutionReaction->rate();
            nDissolutionReactions++;
        }

        localNeighbors += dissolutionReaction->nNeighbors();
    }

    if (localDissolutionRate > 100)
    {
        solver()->exit("Pressure too high?");
    }

    localNeighbors /= m_dissolutionReactions.size();
    if (nDissolutionReactions != 0)
    {
        localDissolutionRate /= nDissolutionReactions;
    }

    m_neighbours += localNeighbors;
    m_dissolutionRate += localDissolutionRate;

    double avgNeighbours = m_neighbours/m_counter;
    double avgDissolutionRate = m_dissolutionRate/m_counter;

    m_eqConc += avgDissolutionRate/(4 - avgNeighbours);

    m_counter++;
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
    double eqConc = m_eqConcEvent->value();

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

    if (m_allConcentrations.size() > 2 && fabs((cNew + m_cPrev)/2 - meanConcentration()) < m_treshold)
    {
        terminateLoop("Concentration converged");
    }

    else
    {
        cout << fabs((cNew + m_cPrev)/2 - meanConcentration())/m_treshold << endl;
        solver()->setTargetConcentration(cNew);
        m_eqConcEvent->restart();

        for (Deposition *depositionReaction : m_depositionReactions)
        {
            depositionReaction->registerUpdateFlag(QuasiDiffusionReaction::UpdateFlags::CALCULATE);
        }

        m_cPrev = cNew;
    }
}
