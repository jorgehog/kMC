#include "boundary.h"

#include "../kmcsolver.h"

#include <climits>

#include <armadillo>



using namespace arma;
using namespace kMC;


Boundary::Boundary(const uint dimension, const uint orientation) :
    m_dimension(dimension),
    m_orientation(orientation)
{

}

Boundary::~Boundary()
{

}


void Boundary::setMainSolver(KMCSolver *solver)
{

    m_solver = solver;

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


uint Boundary::BLOCKED_COORDINATE = (uint)ULLONG_MAX;

KMCSolver* Boundary::m_solver;


