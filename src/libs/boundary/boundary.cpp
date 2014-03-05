#include "boundary.h"

#include "../kmcsolver.h"

#include <armadillo>



using namespace arma;
using namespace kMC;


Boundary::Boundary(const uint dimension, const uint orientation)
{

    switch (dimension) {
    case X:
        m_span = m_NX;
        break;
    case Y:
        m_span = m_NY;
        break;
    case Z:
        m_span = m_NZ;
        break;
    default:
        break;
    }

    m_orientation = orientation;

}

Boundary::~Boundary()
{

}

void Boundary::setMainSolver(KMCSolver *solver)
{

    m_mainSolver = solver;

    m_NX = solver->getNX();
    m_NY = solver->getNY();
    m_NZ = solver->getNZ();

    m_NXYZ = {m_NX, m_NY, m_NZ};

    BLOCKED_COORDINATE = 2*max(uvec{m_NX, m_NY, m_NZ}) + 1;

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

    uvec xyz = {x, y, z};

    for (uint i = 0; i < 3; ++i)
    {
        if (xyz(i) >= m_NXYZ(i))
        {
            loc(i) = 1;
        }

        else
        {
            loc(i) = 0;
        }
    }

}

uint Boundary::BLOCKED_COORDINATE;

uint Boundary::m_NX;
uint Boundary::m_NY;
uint Boundary::m_NZ;

uvec Boundary::m_NXYZ;

KMCSolver* Boundary::m_mainSolver;


