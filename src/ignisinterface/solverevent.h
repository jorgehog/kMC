#pragma once

#include "kmcevent.h"

#include "../kmcsolver.h"

using ignis::Event;

namespace kMC
{

class SolverEvent : public KMCEvent
{
public:

    SolverEvent(KMCSolver *solver) :
        KMCEvent("kmC::SolverEvent"),
        m_solver(solver)
    {

    }

    void initialize()
    {
        m_solver->initialize();
    }

protected:

    void execute()
    {
        m_solver->singleLoop();
    }

private:

    KMCSolver *m_solver;

};


}
