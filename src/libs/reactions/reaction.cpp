#include "reaction.h"
#include "../kmcsolver.h"

uint Reaction::IDcount = 0;

Reaction::Reaction():
    m_ID(IDcount++)
{

}

void Reaction::setMainsolver(KMCSolver *solver)
{
        NX = solver->NX;
        NY = solver->NY;
        NZ = solver->NZ;
        mainSolver = solver;
}
