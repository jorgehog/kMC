#include "reaction.h"

#include "../kmcsolver.h"
#include "site.h"

#include "../debugger/debugger.h"

using namespace kMC;

Reaction::Reaction(Site *currentSite, const string name):
    name(name),
    m_ID(m_IDCount++),
    m_reactionSite(currentSite),
    m_lastUsedEnergy(UNSET_ENERGY),
    m_rate(UNSET_RATE),
    m_updateFlag(UNSET_UPDATE_FLAG)
{

}

Reaction::~Reaction()
{
    m_reactionSite = NULL;
}

const string Reaction::info(int xr, int yr, int zr, string desc) const
{
    stringstream s;
    s << "[" << name << "]:" << "\n";
    s << "   rate: " << m_rate << "  ";
    s << "updateFlags: ";

    for (int flag : m_updateFlags)
    {
        s << flag;
    }

    s << "  ";
    s << "Selected flag: " << m_updateFlag << "  ";
    s << "Blocked? " << !isAllowed() << "\n";
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

    if (!Debugger::enabled) return "";

    int X, Y, Z;
    X = 0;
    Y = 0;
    Z = 0;

    stringstream s;

    const Reaction * lastReaction = Debugger::lastCurrentReaction;

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

void Reaction::setRate(const double rate)
{
    m_lastUsedEnergy = m_reactionSite->energy();
    m_rate = rate;
}

const uint &Reaction::NX()
{
    return m_solver->NX();
}

const uint &Reaction::NY()
{
    return m_solver->NY();
}

const uint &Reaction::NZ()
{
    return m_solver->NZ();
}


void Reaction::setMainSolver(KMCSolver *solver)
{
    m_solver = solver;
}

void Reaction::loadConfig(const Setting &setting)
{

    m_beta            = getSurfaceSetting<double>(setting, "beta");
    m_linearRateScale = getSurfaceSetting<double>(setting, "scale");

}

void Reaction::reset()
{

    m_lastUsedEnergy = UNSET_ENERGY;

    m_rate = UNSET_RATE;

    m_updateFlag = UNSET_UPDATE_FLAG;

}

void Reaction::selectTriumphingUpdateFlag()
{

    m_updateFlag = *std::min_element(m_updateFlags.begin(), m_updateFlags.end());

    KMCDebugger_Assert(m_updateFlag, !=, UNSET_UPDATE_FLAG, "Update flag was not initialized correctly.", info());

    m_updateFlags.clear();

}

const double   Reaction::UNSET_RATE = -1337;
const double   Reaction::UNSET_ENERGY = -13371337;

KMCSolver*     Reaction::m_solver;

double         Reaction::m_beta;
double         Reaction::m_linearRateScale;

uint           Reaction::m_IDCount = 0;


ostream & operator << (ostream& os, const Reaction& ss)
{
    os << ss.str();
    return os;
}
