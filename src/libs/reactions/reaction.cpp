#include "reaction.h"

#include "../kmcsolver.h"
#include "site.h"

#include "../debugger/kmcdebugger.h"

Reaction::Reaction(string name):
    name(name),
    m_ID(IDcount++),
    m_rate(UNSET_RATE),
    m_updateFlag(UNSET_UPDATE_FLAG)
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

const uint &Reaction::x() const
{
    return m_reactionSite->x();
}

const uint &Reaction::y() const
{
    return m_reactionSite->y();
}

const uint &Reaction::z() const
{
    return m_reactionSite->z();
}


string Reaction::getFinalizingDebugMessage() const
{
#ifndef KMC_NO_DEBUG

    if (!KMCDebugger::enabled) return "";

    int X, Y, Z;
    X = 0;
    Y = 0;
    Z = 0;

    stringstream s;

    const Reaction * lastReaction = KMCDebugger::lastCurrentReaction;

    if (lastReaction != NULL)
    {
        m_reactionSite->distanceTo(lastReaction->reactionSite(), X, Y, Z);
    }
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

void Reaction::selectTriumphingUpdateFlag()
{

    m_updateFlag = *std::min_element(m_updateFlags.begin(), m_updateFlags.end());

    KMCDebugger_Assert(m_updateFlag, !=, UNSET_UPDATE_FLAG, "Update flag was not initialized correctly.", info());

    m_updateFlags.clear();

}

const double   Reaction::UNSET_RATE = -1.337;

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
