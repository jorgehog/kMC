#include "boundary.h"

#include "../kmcsolver.h"
#include "../site.h"

#include "../debugger/debugger.h"

#include <limits>

#include <armadillo>



using namespace arma;
using namespace kMC;


Boundary::Boundary(const uint dimension, const uint orientation, const uint type) :
    type(type),
    m_dimension(dimension),
    m_orientation(orientation)
{



}

Boundary::~Boundary()
{
    m_boundarySites.clear();
}

void Boundary::setupBoundarySites()
{

    uint xi = (orientation() == 0) ? 0 : (span() - 1);

    m_boundarySites.clear();

    if (dimension() == X)
    {

        for (uint y = 0; y < NY(); ++y)
        {
            for (uint z = 0; z < NZ(); ++z)
            {
                m_boundarySites.push_back(solver()->getSite(xi, y, z));
            }
        }

    }

    else if (dimension() == Y)
    {

        for (uint x = 0; x < NX(); ++x)
        {
            for (uint z = 0; z < NZ(); ++z)
            {
                m_boundarySites.push_back(solver()->getSite(x, xi, z));
            }
        }

    }

    else
    {

        KMCDebugger_Assert(dimension(), ==, Z, "This else should always correspond to dim=2");

        for (uint x = 0; x < NX(); ++x)
        {
            for (uint y = 0; y < NY(); ++y)
            {
                m_boundarySites.push_back(solver()->getSite(x, y, xi));
            }
        }

    }

    KMCDebugger_Assert(m_boundarySites.size(), ==, NX()*NY()*NZ()/span(), "mismatch in boundary site setup.");

}

void Boundary::distanceFromSite(const Site *site, int &dxi, bool abs)
{
    uint xi = orientation() == 0 ? 0 : span() - 1;

    dxi = getDistanceBetween(site->r(dimension()), xi);

    if (abs)
    {
        dxi = std::abs(dxi);
    }

}


void Boundary::setMainSolver(KMCSolver *solver)
{
    m_solver = solver;
    m_currentBoundaries.resize(3);
}

bool Boundary::isCompatible(const int type1, const int type2, bool reverse)
{
    bool compatible = !(type1 == Periodic && type2 != Periodic);

    if (reverse)
    {
        compatible = (compatible && isCompatible(type2, type1, false));
    }

    return compatible;

}

void Boundary::setupLocations(const uint x, const uint y, const uint z, uvec3 &loc)
{
    //make x, y, z boundary static site members?

    uvec xyz = {x, y, z};

    for (uint i = 0; i < 3; ++i)
    {
        if (xyz(i) >= N(i)/2)
        {
            loc(i) = 1;
        }

        else
        {
            loc(i) = 0;
        }
    }

}


const uint & Boundary::NX()
{
    return m_solver->NX();
}

const uint & Boundary::NY()
{
    return m_solver->NY();
}

const uint & Boundary::NZ()
{
    return m_solver->NZ();
}

uint Boundary::N(const uint i)
{
    return m_solver->N(i);
}

void Boundary::setupCurrentBoundaries(const uint x, const uint y, const uint z)
{
    uvec3 loc;
    setupLocations(x, y, z, loc);

    m_currentBoundaries.at(0) = Site::boundaries(0, loc(0));
    m_currentBoundaries.at(1) = Site::boundaries(1, loc(1));
    m_currentBoundaries.at(2) = Site::boundaries(2, loc(2));

}


uint Boundary::BLOCKED_COORDINATE = std::numeric_limits<uint>::max();

KMCSolver* Boundary::m_solver;


vector<const Boundary*> Boundary::m_currentBoundaries;
