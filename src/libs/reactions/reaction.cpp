#include "reaction.h"
#include "../kmcsolver.h"


Reaction::Reaction():
    m_ID(IDcount++)
{

}

Reaction::~Reaction()
{

}


void Reaction::setSolverPtr(KMCSolver *solver)
{
    NX = solver->getNX();
    NY = solver->getNY();
    NZ = solver->getNZ();

    mainSolver = solver;
}

void Reaction::loadReactionSettings(const Setting &setting)
{
    beta = getSurfaceSetting<double>(setting, "beta");
    mu   = getSurfaceSetting<double>(setting, "scale");
}

KMCSolver *Reaction::mainSolver;
double Reaction::beta;
double Reaction::mu;

uint Reaction::NX;
uint Reaction::NY;
uint Reaction::NZ;

uint Reaction::IDcount = 0;
