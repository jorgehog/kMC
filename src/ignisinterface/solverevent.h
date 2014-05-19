#pragma once

#include "kmcevent.h"

#include "../kmcsolver.h"

#include "../soluteparticle.h"
#include "../potential/potential.h"

using ignis::Event;

namespace kMC
{

class SolverEvent : public KMCEvent
{
public:

    SolverEvent() : KMCEvent("kmC::SolverEvent")
    {

    }

    void initialize()
    {
        KMCDebugger_Init();

        KMCDebugger_AssertBool(Site::boundariesIsInitialized(), "Boundaries are not initialized.");
        KMCDebugger_Assert(Site::_refCount(), !=, 0, "Sites need to be initialized.");

        if (SoluteParticle::ss != NULL)
        {
            SoluteParticle::ss->initialize();
        }

        solver()->getRateVariables();

        m_totalTime = 0;

        resetReaction();

    }

    const double & totalTime() const
    {
        return m_totalTime;
    }

    const Reaction * selectedReaction() const
    {
        return m_selectedReaction;
    }

    void resetReaction()
    {
        m_selectedReaction = NULL;
    }

protected:

    void execute()
    {
        R = solver()->kTot()*KMC_RNG_UNIFORM();

        choice = solver()->getReactionChoice(R);

        m_selectedReaction = solver()->allPossibleReactions().at(choice);
        KMCDebugger_SetActiveReaction(m_selectedReaction);

        m_selectedReaction->execute();

        Site::updateBoundaries();

        solver()->getRateVariables();

        m_totalTime += Reaction::linearRateScale()/solver()->kTot();

        //To counter buildup of roundoff errors
        if (nTimesExecuted % 10000 == 0)
        {
            solver()->remakeAccuAllRates();
        }

    }

private:

    double R;

    double m_totalTime;

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

protected:

    virtual void dumpFile() const = 0;

    void execute()
    {
        if (nTimesExecuted%MainLattice::nCyclesPerOutput != 0)
        {
            return;
        }

        cout << "Storing file: " << m_outputCounter << endl;
        dumpFile();

        m_outputCounter++;
    }

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
