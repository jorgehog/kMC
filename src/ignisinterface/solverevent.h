#pragma once

#include "kmcevent.h"

#include "../kmcsolver.h"

using ignis::Event;

namespace kMC
{

class SolverEvent : public KMCEvent
{
public:

    SolverEvent() :
        KMCEvent("kmC::SolverEvent")
    {

    }

    void initialize()
    {
        KMCDebugger_Init();

        KMCDebugger_AssertBool(Site::boundariesIsInitialized(), "Boundaries are not initialized.");
        KMCDebugger_Assert(Site::_refCount(), !=, 0, "Sites need to be initialized.");

        solver()->getRateVariables();

        totalTime = 0;

    }

protected:

    void execute()
    {
        R = solver()->kTot()*KMC_RNG_UNIFORM();

        choice = solver()->getReactionChoice(R);

        selectedReaction = solver()->allPossibleReactions().at(choice);
        KMCDebugger_SetActiveReaction(selectedReaction);

        selectedReaction->execute();

        Site::updateBoundaries();

        solver()->getRateVariables();

        totalTime += Reaction::linearRateScale()/solver()->kTot();

    }

private:

    double R;

    double totalTime;

    uint choice;

    Reaction * selectedReaction;

};


class DumpXYZ : public KMCEvent
{
public:

    DumpXYZ() : KMCEvent("DumpXYZ") {}

protected:

    void initialize()
    {
        outputCounter = 0;
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

};

}
