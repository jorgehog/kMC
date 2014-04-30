#pragma once

#include "kmcevent.h"

#include "../kmcsolver.h"

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

    }

private:

    double R;

    double m_totalTime;

    uint choice;

    Reaction * m_selectedReaction;

};


class DumpXYZ : public KMCEvent
{
public:

    DumpXYZ() : KMCEvent("DumpXYZ"), offset(0) {}

    void setOffset(const uint offset)
    {
        this->offset = offset;
    }

protected:

    void initialize()
    {
        outputCounter = offset;
    }

    void execute()
    {
        if (nTimesExecuted%MainLattice::nCyclesPerOutput != 0)
        {
            return;
        }

        cout << "Storing XYZ: " << outputCounter << endl;
        solver()->dumpXYZ(outputCounter++);

    }

private:

    uint outputCounter;
    uint offset;

};

}
