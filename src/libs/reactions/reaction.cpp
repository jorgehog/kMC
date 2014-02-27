#include "reaction.h"
#include "../kmcsolver.h"

#include "../debugger/kmcdebugger.h"

Reaction::Reaction(string name):
    name(name),
    m_ID(IDcount++),
    m_rate(UNSET_RATE),
    m_updateFlag(defaultUpdateFlag)
{

}

Reaction::~Reaction()
{

}

const string Reaction::info(int xr, int yr, int zr, string desc) const
{
    stringstream s;
    s << "[" << name << " " << m_ID << "/" << IDcount << "]:" << "\n";
    s << "   rate: " << m_rate << "  ";
    s << "updateFlags: ";

    for (int flag : m_updateFlags)
    {
        s << flag;
    }

    s << "  ";
    s << "Selected flag: " << m_updateFlag << "  ";
    s << "Blocked? " << !isNotBlocked() << "\n";
    s << "@";
    s << m_reactionSite->info(xr, yr, zr, desc);
    s << "\n";

    return s.str();

}

string Reaction::getFinalizingDebugMessage() const
{
#ifndef KMC_NO_DEBUG
    int X, Y, Z;
    stringstream s;

    const Reaction * lastReaction = KMCDebugger_GetReaction(lastCurrent);
    const Site* site = lastReaction->reactionSite();

    m_reactionSite->distanceTo(site, X, Y, Z);

    s << info();
    s << "\nLast active reaction site marked on current site:\n\n";
    s << m_reactionSite->info(X, Y, Z);

    return s.str();
#else
    return "";
#endif
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

void Reaction::getTriumphingUpdateFlag()
{
    if (m_updateFlags.empty())
    {
        m_updateFlag = defaultUpdateFlag;
    }

    else
    {
        m_updateFlag = *std::min_element(m_updateFlags.begin(), m_updateFlags.end());
    }

    clearUpdateFlags();

}

const double   Reaction::UNSET_RATE = -1;

KMCSolver*     Reaction::mainSolver;

double         Reaction::beta;
double         Reaction::m_linearRateScale;

uint           Reaction::NX;
uint           Reaction::NY;
uint           Reaction::NZ;

uint           Reaction::IDcount = 0;


ostream & operator << (ostream& os, const Reaction& ss)
{
    os << ss.str();
    return os;
}
