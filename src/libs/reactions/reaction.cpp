#include "reaction.h"
#include "../kmcsolver.h"


Reaction::Reaction(string name):
    name(name),
    m_ID(IDcount++),
    m_rate(UNSET_RATE)
{

}

Reaction::~Reaction()
{

}

void Reaction::dumpInfo(int xr, int yr, int zr) const
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
    os << ss.name << "@(" << ss.x() << "," << ss.y() << "," << ss.z() << ") [" << ss.getInfoSnippet() << "]";
    return os;
}
