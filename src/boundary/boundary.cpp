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
    m_orientation(orientation),
    m_bound((orientation == 0) ? 0 : (span() - 1)),
    m_initialized(false)
{



}

Boundary::~Boundary()
{

}

uint Boundary::transformCoordinate(const int xi) const
{
    if ((xi >= (int)span()) || (xi < 0))
    {
        return BLOCKED_COORDINATE;
    }

    else
    {
        return xi;
    }

}

void Boundary::getBoundarySite(uint n, uint &x, uint&y, uint &z)
{
    if (dimension() == Z)
    {
        x = NX()%n;
        y = n - x;
        z = m_bound;

    }
    else if (dimension() == Y)
    {
        x = NX()%n;
        y = m_bound;
        z = n - x;
    }
    else
    {
        KMCDebugger_Assert(dimension(), ==, X);

        x = m_bound;
        y = NY()%n;
        z = n - x;
    }

}

void Boundary::distanceFrom(const uint xi, int &dxi, bool abs)
{
    uint loc = orientation() == 0 ? 0 : span() - 1;

    dxi = getDistanceBetween(xi, loc);

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

uint Boundary::getLocation(const uint xi, const uint dim, const uint shift)
{
    if (xi >= N(dim)/2 + shift)
    {
        return 1;
    }

    return 0;

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

void Boundary::setupCurrentBoundary(const uint x, const uint dim, const uint shift)
{
    m_currentBoundaries.at(dim) = Site::boundaries(dim, getLocation(x, dim, shift));
}

void Boundary::setupCurrentBoundaries(const uint x, const uint y, const uint z, const uint shift)
{
    setupCurrentBoundary(x, 0, shift);
    setupCurrentBoundary(y, 1, shift);
    setupCurrentBoundary(z, 2, shift);
}


uint Boundary::BLOCKED_COORDINATE = std::numeric_limits<uint>::max();

KMCSolver* Boundary::m_solver;


vector<const Boundary*> Boundary::m_currentBoundaries;
