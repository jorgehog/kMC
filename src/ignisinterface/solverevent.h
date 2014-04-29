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
