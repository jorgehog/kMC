#include "reaction.h"
#include "../kmcsolver.h"


Reaction::Reaction():
    m_ID(IDcount++)
{

}

void Reaction::setSolverPtr(KMCSolver *solver)
{
    NX = solver->getNX();
    NY = solver->getNY();
    NZ = solver->getNZ();

    mainSolver = solver;
}

void Reaction::loadTemperature(const Setting &setting)
{
    beta = getSurfaceSetting<double>(setting, "beta");
}

KMCSolver *Reaction::mainSolver;
double Reaction::beta;

uint Reaction::NX;
uint Reaction::NY;
uint Reaction::NZ;

uint Reaction::IDcount = 0;
