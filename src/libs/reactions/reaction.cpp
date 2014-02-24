#include "reaction.h"
#include "../kmcsolver.h"

#include "../debugger/kmcdebugger.h"

Reaction::Reaction(string name):
    name(name),
    m_ID(IDcount++),
<<<<<<< HEAD
    m_siteReactionArrayIndex(UNSET_ARRAY_INDEX),
    m_rate(UNSET_RATE)
=======
    m_rate(UNSET_RATE),
    m_updateFlag(UNSET_UPDATE_FLAG)
>>>>>>> experimental2
{

}

Reaction::~Reaction()
{

}

const string Reaction::info(int xr, int yr, int zr, string desc) const
{
    stringstream s;
    s << "[" << name << " " << m_ID << "/" << IDcount << "]:" << "\n";
    s << "rate: " << m_rate << "\n";
    s << "updateFlags: ";

    for (int flag : m_updateFlags)
    {
        s << flag;
    }

    s << "\n";
    s << "Selected flag: " << m_updateFlag << "\n";

    s << "@{" << "\n";
    s << m_reactionSite->info(xr, yr, zr, desc);
    s << "\n}";

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
    s << "\nLast active reaction site marked on current site:\n";
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

<<<<<<< HEAD
void Reaction::initialize()
{
    if (m_reactionSite->isActive() && isNotBlocked())
    {
        assert(false && "NOT TESTED");
        enable();
        calcRate();
    }
}

void Reaction::update()
{

    if (!m_reactionSite->isActive())
    {
        m_siteReactionArrayIndex = UNSET_ARRAY_INDEX;

        return;
    }

    if (isNotBlocked())
    {
        enable();
=======
void Reaction::getTriumphingUpdateFlag()
{
    if (m_updateFlags.empty())
    {
        m_updateFlag = defaultUpdateFlag;
>>>>>>> experimental2
    }

    else
    {
<<<<<<< HEAD
        disable();
    }

    calcRate();

}
void Reaction::enable()
{
    m_reactionSite->enableReaction(this);
}

void Reaction::disable()
{
    m_reactionSite->disableReaction(this);
}

const double Reaction::UNSET_RATE = -1;
const uint   Reaction::UNSET_ARRAY_INDEX = 27;
=======
        m_updateFlag = *std::min_element(m_updateFlags.begin(), m_updateFlags.end());
    }

    clearUpdateFlags();

}

const double   Reaction::UNSET_RATE = -1;
>>>>>>> experimental2

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
