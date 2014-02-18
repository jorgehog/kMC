#include "reaction.h"
#include "../kmcsolver.h"


Reaction::Reaction():
    m_ID(IDcount++),
    m_rate(UNSET_RATE)
{

}

Reaction::~Reaction()
{

}

void Reaction::dumpInfo(int xr, int yr, int zr)
{

    cout << "[Reaction " << m_ID << "/" << IDcount << "]:" << endl;
    cout << "rate: " << m_rate << endl;
    cout << "@{" << endl;
    m_reactionSite->dumpInfo(xr, yr, zr);
    cout << "\n}" << endl;

}


void Reaction::setSolverPtr(KMCSolver *solver)
{

    NX = solver->getNX();
    NY = solver->getNY();
    NZ = solver->getNZ();

    mainSolver = solver;

}

void Reaction::loadConfig(const Setting &setting)
{

    beta              = getSurfaceSetting<double>(setting, "beta");
    m_linearRateScale = getSurfaceSetting<double>(setting, "scale");

}

const double Reaction::UNSET_RATE = -1;

KMCSolver*   Reaction::mainSolver;

double       Reaction::beta;
double       Reaction::m_linearRateScale;

uint         Reaction::NX;
uint         Reaction::NY;
uint         Reaction::NZ;

uint         Reaction::IDcount = 0;
