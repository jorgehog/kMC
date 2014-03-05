#include "boundary.h"

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

}


uint Boundary::m_NX;
uint Boundary::m_NY;
uint Boundary::m_NZ;

