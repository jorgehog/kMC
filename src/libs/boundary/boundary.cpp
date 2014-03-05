#include "boundary.h"

#include <armadillo>

using namespace arma;
using namespace kMC;

Boundary::Boundary()
{
}

void Boundary::loadConfig(const Setting &setting)
{

    const Setting & BoxSize = getSurfaceSetting(setting, "BoxSize");


    m_NX = BoxSize[0];
    m_NY = BoxSize[1];
    m_NZ = BoxSize[2];

    BLOCKED_COORDINATE = max(uvec{m_NX, m_NY, m_NZ}) + 1;

}

uint Boundary::BLOCKED_COORDINATE;

uint Boundary::m_NX;
uint Boundary::m_NY;
uint Boundary::m_NZ;


