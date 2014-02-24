#include "reaction.h"
#include "../kmcsolver.h"


Reaction::Reaction(string name):
    name(name),
    m_ID(IDcount++),
    m_rate(UNSET_RATE),
    m_type(std)
{

}

Reaction::~Reaction()
{

}

const string Reaction::info(int xr, int yr, int zr, string desc) const
{
    stringstream s;
    s << "[" << name << " " << m_ID << "/" << IDcount << "]:" << endl;
    s << "rate: " << m_rate << endl;
    s << "@{" << endl;
    s << m_reactionSite->info(xr, yr, zr, desc);
    s << "\n}" << endl;

    return s.str();

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

const double Reaction::UNSET_RATE = -1;

KMCSolver*   Reaction::mainSolver;

double       Reaction::beta;
double       Reaction::m_linearRateScale;

uint         Reaction::NX;
uint         Reaction::NY;
uint         Reaction::NZ;

uint         Reaction::IDcount = 0;


ostream & operator << (ostream& os, const Reaction& ss)
{
    os << ss.str();
    return os;
}
