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

        stringstream s;
        s << "kMC" << outputCounter++ << ".xyz";

        ofstream o;
        o.open("outfiles/" + s.str());



        stringstream surface;
        stringstream crystal;
        stringstream solution;

        s.str(string());

        for (SoluteParticle *particle : solver()->particles())
        {

            s << "\n"
              << particle->particleStateShortName() << " "
              << particle->x() << " " << particle->y() << " " << particle->z() << " "
              << particle->nNeighborsSum() << " "
              << particle->energy();

            if (particle->isSurface())
            {
                surface << s.str();
            }

            else if (particle->isCrystal())
            {
                crystal << s.str();
            }

            else
            {
                solution << s.str();
            }

            s.str(string());

        }

        o << solver()->particles().size() << "\n - " << surface.str() << crystal.str() << solution.str();
        o.close();

    }

private:

    uint outputCounter;

};

}
