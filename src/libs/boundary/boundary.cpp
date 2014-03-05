#include "boundary.h"

#include "../kmcsolver.h"

#include <armadillo>

using namespace arma;
using namespace kMC;

Boundary::Boundary()
{

}

void Boundary::setMainSolver(KMCSolver *solver)
{

    m_mainSolver = solver;

    m_NX = solver->getNX();
    m_NY = solver->getNY();
    m_NZ = solver->getNZ();

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

int  Boundary::BLOCKED_COORDINATE;

uint Boundary::m_NX;
uint Boundary::m_NY;
uint Boundary::m_NZ;

KMCSolver* Boundary::m_mainSolver;


