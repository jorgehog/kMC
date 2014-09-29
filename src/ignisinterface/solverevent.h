#pragma once

#include "kmcevent.h"

#include "../kmcsolver.h"

#include "../soluteparticle.h"
#include "../potential/potential.h"

#include "../debugger/debugger.h"

#include "../reactions/reaction.h"

#include <BADAss/badass.h>

using ignis::Event;

namespace kMC
{

class SolverEvent : public KMCEvent
{
public:

    SolverEvent() : KMCEvent("kMC::SolverEvent")
    {

    }

    void initialize()
    {
        KMCDebugger_Init();

        BADAssBool(Site::boundariesIsInitialized(), "Boundaries are not initialized.");
        BADAss(Site::_refCount(), !=, 0, "Sites need to be initialized.");

        if (SoluteParticle::ss != NULL)
        {
            SoluteParticle::ss->initialize();
        }

        solver()->getRateVariables();

        m_totalTime = 0;

        resetReaction();

        if (!solver()->localUpdating())
        {
            solver()->remakeAccuAllRates();
        }

    }

    const double & totalTime() const
    {
        return m_totalTime;
    }

    const double &lastTimeStep() const
    {
        return m_lastTimeStep;
    }

    const Reaction * selectedReaction() const
    {
        return m_selectedReaction;
    }

    void resetReaction()
    {
        m_selectedReaction = NULL;
    }

    void reset()
    {
        Site::updateBoundaries();

        solver()->getRateVariables();

        //To counter buildup of roundoff errors
        if ((m_cycle % 10000 == 0) && solver()->localUpdating())
        {
            solver()->remakeAccuAllRates();
        }

        BADAssBool(!solver()->allPossibleReactions().empty(), "No available reactions.");
        BADAss(solver()->allPossibleReactions().size(), ==, solver()->accuAllRates().size(), "These vectors should be equal of length.");
        BADAssClose(solver()->accuAllRates().at(0), solver()->allPossibleReactions().at(0)->rate(), solver()->minRateThreshold(), "zeroth accuallrate should be the first rate.");
        BADAssClose(solver()->accuAllRates().at(solver()->accuAllRates().size()- 1), solver()->kTot(), solver()->minRateThreshold(), "kTot should be the last element of accuAllRates");

    }

    void execute()
    {
        R = solver()->kTot()*KMC_RNG_UNIFORM();

        choice = solver()->getReactionChoice(R);

        m_selectedReaction = solver()->allPossibleReactions().at(choice);
        KMCDebugger_SetActiveReaction(m_selectedReaction);

        m_selectedReaction->execute();

        m_lastTimeStep = -Reaction::linearRateScale()*std::log(KMC_RNG_UNIFORM())/solver()->kTot();

        BADAss(m_lastTimeStep, >, 0, "timestep should be posiive.");

        m_totalTime += m_lastTimeStep;
    }

private:

    double R;

    double m_totalTime;

    double m_lastTimeStep;

    uint choice;

    Reaction * m_selectedReaction;

};

class DumpFile : public KMCEvent
{

public:

    DumpFile() : KMCEvent("DumpFile"), m_offset(0) {}

    void setOffset(const uint offset)
    {
        m_offset = offset;
    }

    void initialize()
    {
        m_outputCounter = m_offset;
    }

    void execute()
    {
        if (m_cycle%solver()->mainLattice()->outputSpacing() != 0)
        {
            return;
        }

        cout << "Storing file: " << m_outputCounter << endl;
        dumpFile();

        m_outputCounter++;
    }

protected:

    virtual void dumpFile() const = 0;

    const uint &outputCounter() const
    {
        return m_outputCounter;
    }


private:

    uint m_outputCounter;
    uint m_offset;

};

class DumpXYZ : public DumpFile
{
public:

    DumpXYZ() : DumpFile() {}

protected:

    void dumpFile() const
    {
        solver()->dumpXYZ(outputCounter());
    }

};


class DumpLAMMPS: public DumpFile
{
public:

    DumpLAMMPS() : DumpFile() {}

protected:

    void dumpFile() const
    {
        solver()->dumpLAMMPS(outputCounter());
    }

};

}
