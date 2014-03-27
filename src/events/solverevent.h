#pragma once

#include "../kmcsolver.h"
#include <ignis.h>

using ignis::Event;

namespace kMC
{


class SolverEvent : public Event
{
public:

    SolverEvent(KMCSolver *solver) :
        Event("kmC::SolverEvent"),
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
