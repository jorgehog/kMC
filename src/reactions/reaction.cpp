#include "reaction.h"

#include "diffusion/diffusionreaction.h"

#include "../kmcsolver.h"
#include "site.h"

#include "../debugger/debugger.h"

using namespace kMC;

Reaction::Reaction(Site *currentSite):
    m_reactionSite(currentSite),
    m_lastUsedEnergy(UNSET_ENERGY),
    m_rate(UNSET_RATE),
    m_updateFlag(UNSET_UPDATE_FLAG),
    m_address(UNSET_ADDRESS)
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
    s << propertyString() << "\n";
    s << "@";
    s << m_reactionSite->info(xr, yr, zr, desc);
    s << "\n";

    return s.str();

}

void Reaction::setLastUsedEnergy()
{
    m_lastUsedEnergy = m_reactionSite->energy();
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

bool Reaction::hasVacantStatus() const
{
    return solver()->isEmptyAddress(m_address);
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

string Reaction::propertyString() const
{
    stringstream s;

    s << "Rate: "          << unsetIf(m_rate,       UNSET_RATE)        << "  ";
    s << "Address: "       << unsetIf(m_address,    UNSET_ADDRESS)     << "  ";
    s << "Selected flag: " << unsetIf(m_updateFlag, UNSET_UPDATE_FLAG) << "  ";

    s << "Allowed? " << !isAllowed() << "  ";
    s << "Vacant? "  << hasVacantStatus();

    return s.str();
}

void Reaction::setRate(const double rate)
{

    KMCDebugger_Assert(rate, !=, 0, "This reaction should be deactive.", getFinalizingDebugMessage());

    m_solver->registerReactionChange(this, rate);

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

void Reaction::setBeta(const double beta)
{

    DiffusionReaction::setBetaChangeScaleFactor(std::exp(beta - m_beta));

    for (Reaction * r : m_solver->allPossibleReactions())
    {
        r->registerBetaChange(beta);
    }

    solver()->onAllRatesChanged();

    m_beta = beta;

}

void Reaction::reset()
{

    m_lastUsedEnergy = UNSET_ENERGY;

    m_rate = UNSET_RATE;

    m_updateFlag = UNSET_UPDATE_FLAG;

    m_address = UNSET_ADDRESS;

}

bool Reaction::isAllowedAndActive() const
{
    return isAllowed() && reactionSite()->isActive();
}

const string Reaction::name = "Reaction";

KMCSolver*   Reaction::m_solver;

double       Reaction::m_beta = 1.0;
double       Reaction::m_linearRateScale = 1.0;

uint         Reaction::m_IDCount = 0;

const double Reaction::UNSET_RATE     = -1337;

const double Reaction::UNSET_ENERGY   = -13371337;

const uint   Reaction::UNSET_ADDRESS  = std::numeric_limits<uint>::max();


ostream & operator << (ostream& os, const Reaction& ss)
{
    os << ss.str();
    return os;
}
