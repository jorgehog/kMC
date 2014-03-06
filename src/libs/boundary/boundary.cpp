#include "boundary.h"

#include "../kmcsolver.h"

#include "../debugger/debugger.h"

#include <armadillo>



using namespace arma;
using namespace kMC;


Boundary::Boundary(const uint dimension, const uint orientation) :
    m_span(m_NXYZ(dimension)),
    m_dimension(dimension),
    m_orientation(orientation)
{

}

Boundary::~Boundary()
{

}

void Boundary::initialize()
{
    uint xi = 0;

    if (orientation() == 1)
    {
        xi = span() - 1;
    }

    if (dimension() == X)
    {

        for (uint y = 0; y < NY(); ++y) {
            for (uint z = 0; z < NZ(); ++z) {
                if (!mainSolver()->getSite(xi, y, z)->isActive())
                {

                    mainSolver()->getSite(xi, y, z)->spawnAsFixedCrystal();
                }
            }
        }

    }

    else if (dimension() == Y)
    {

        for (uint x = 0; x < NX(); ++x) {
            for (uint z = 0; z < NZ(); ++z) {
                if (!mainSolver()->getSite(x, xi, z)->isActive())
                {
                    mainSolver()->getSite(x, xi, z)->spawnAsFixedCrystal();
                }
            }
        }

    }

    else
    {
        KMCDebugger_Assert(dimension(), ==, Z, "This else should always correspond to dim=2");

        for (uint x = 0; x < NX(); ++x) {
            for (uint y = 0; y < NY(); ++y) {
                if (!mainSolver()->getSite(x, y, xi)->isActive())
                {
                    mainSolver()->getSite(x, y, xi)->spawnAsFixedCrystal();
                }
            }
        }

    }


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


