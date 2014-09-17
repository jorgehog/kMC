#include "reaction.h"

#include "diffusion/diffusionreaction.h"

#include "../kmcsolver.h"
#include "../soluteparticle.h"

#include "../ignisinterface/solverevent.h"

#include <BADAss/badass.h>

using namespace kMC;

Reaction::Reaction(SoluteParticle *reactant):
    m_reactant(reactant),
    m_lastUsedEnergy(UNSET_ENERGY),
    m_rate(UNSET_RATE),
    m_updateFlag(UNSET_UPDATE_FLAG),
    m_address(UNSET_ADDRESS)
{
    refCount++;
}

Reaction::~Reaction()
{

    disable();

    m_reactant = NULL;

    BADAss(refCount, !=, 0);

    if (solver()->solverEvent()->selectedReaction() == this)
    {
        solver()->resetLastReaction();
    }

    refCount--;
}

const string Reaction::info(int xr, int yr, int zr, string desc) const
{
    stringstream s;
    s << "[" << name() << "]:" << "\n";
    s << propertyString() << "\n";
    s << "@";
    s << m_reactant->info(xr, yr, zr, desc);
    s << "\n";

    return s.str();

}

void Reaction::setLastUsedEnergy()
{
    m_lastUsedEnergy = m_reactant->energy();
}


const uint &Reaction::x() const
{
    return m_reactant->x();
}

const uint &Reaction::y() const
{
    return m_reactant->y();
}

const uint &Reaction::z() const
{
    return m_reactant->z();
}

bool Reaction::hasVacantStatus() const
{
    return solver()->isEmptyAddress(m_address);
}

string Reaction::propertyString() const
{
    stringstream s;

    s << "Rate: "          << unsetIf(m_rate,       UNSET_RATE)        << "  ";
    s << "Address: "       << unsetIf(m_address,    UNSET_ADDRESS)     << "  ";
    s << "Selected flag: " << unsetIf(m_updateFlag, UNSET_UPDATE_FLAG) << "  ";

    s << "Allowed? " <<  isAllowed() << "  ";
    s << "Vacant? "  << hasVacantStatus();

    return s.str();
}

void Reaction::setRate(const double rate)
{

    BADAss(rate, !=, 0, "This reaction should be deactive.", KMCBAI(info()));

    solver()->registerReactionChange(this, rate);

    m_rate = rate;

}


const uint &Reaction::NX()
{
    return solver()->NX();
}

const uint &Reaction::NY()
{
    return solver()->NY();
}

const uint &Reaction::NZ()
{
    return solver()->NZ();
}


void Reaction::setMainSolver(KMCSolver *solver)
{
    m_solver = solver;
}

void Reaction::loadConfig(const Setting &setting)
{

    m_beta            = getSetting<double>(setting, "beta");
    m_linearRateScale = getSetting<double>(setting, "scale");

}

void Reaction::setBeta(const double beta)
{

    if (beta == m_beta)
    {
        return;
    }

    DiffusionReaction::setBetaChangeScaleFactor(std::exp(beta - m_beta));

    for (Reaction * r : solver()->allPossibleReactions())
    {
        r->registerBetaChange(beta);
    }

    solver()->remakeAccuAllRates();

    m_beta = beta;

}

void Reaction::reset()
{

    m_lastUsedEnergy = UNSET_ENERGY;

    m_rate = UNSET_RATE;

    m_updateFlag = UNSET_UPDATE_FLAG;

    m_address = UNSET_ADDRESS;

}

KMCSolver*   Reaction::m_solver;

double       Reaction::m_beta = 1.0;
double       Reaction::m_linearRateScale = 1.0;

uint         Reaction::refCount = 0;

const double Reaction::UNSET_RATE     = -1337;

const double Reaction::UNSET_ENERGY   = -13371337;

const uint   Reaction::UNSET_ADDRESS  = std::numeric_limits<uint>::max();


ostream & operator << (ostream& os, const Reaction& ss)
{
    os << ss.str();
    return os;
}
