#include "boundary.h"

#include "../kmcsolver.h"
#include "../site.h"
#include "../soluteparticle.h"

#include <limits>

#include <armadillo>
#include <BADAss/badass.h>


using namespace arma;
using namespace kMC;


Boundary::Boundary(const uint dimension, const uint orientation, const uint type) :
    m_type(type),
    m_dimension(dimension),
    m_orientation(orientation),
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

void Boundary::getBoundarySite(uint n, uint &x, uint&y, uint &z) const
{
    if (dimension() == Z)
    {
        x = n/NY();
        y = n - x*NY();
        z = interfaceValue(x, y);

    }
    else if (dimension() == Y)
    {
        x = n/NZ();
        z = n - x*NZ();
        y = interfaceValue(x, z);
    }
    else
    {
        BADAss(dimension(), ==, X);

        y = n/NZ();
        z = n - y*NZ();
        x = interfaceValue(y, z);
    }

}

int Boundary::distanceFrom(const uint xi, bool abs) const
{
    int dxi = getDistanceBetween(xi, bigbound());

    if (abs)
    {
        return std::abs(dxi);
    }

    return dxi;

}


void Boundary::setMainSolver(KMCSolver *solver)
{
    m_solver = solver;
    m_currentBoundaries.resize(3);
}

bool Boundary::isCompatible(const int type1, const int type2, bool reverse)
{
    bool compatible = !(type1 == Periodic && type2 != Periodic);
    compatible = compatible && !(type1 == SphericalEdge && type2 != SphericalEdge);

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

void Boundary::performConcentrationBoundaryConditionStep()
{
    counter++;
    if (counter%m_coolDown != 0)
    {
        return;
    }

    BADAss(m_maxEventsPrCycle, <=, boundarySize(), "Max events pr cycle cannot exceed the number of boundary sites.");


    Site * currentSite;

    uint x, y, z;
    uint c = 0;
    uint ce = 0;

    uint targetN = solver()->targetConcentration()*SoluteParticle::getCurrentSolvantVolume();

    if (targetN == 0)
    {
        targetN = 1;
    }

    if (SoluteParticle::nSolutionParticles() > targetN)
    {

        while (SoluteParticle::nSolutionParticles() != targetN && c != boundarySize() && ce != m_maxEventsPrCycle)
        {

            getBoundarySite(KMC_RNG_UNIFORM()*boundarySize(), x, y, z);

            currentSite = solver()->getSite(x, y, z);

            if (currentSite != NULL)
            {
                if (currentSite->isActive())
                {
                    solver()->despawnParticle(currentSite->associatedParticle());
                    ce++;
                }
            }

            c++; //*giggle*
        }

    }

    else
    {

        bool spawned;

        while (SoluteParticle::nSolutionParticles() != targetN && ce != m_maxEventsPrCycle)
        {

            spawned = false;

            SoluteParticle *particle = new SoluteParticle();

            while (c != boundarySize() && !spawned)
            {

                getBoundarySite(KMC_RNG_UNIFORM()*boundarySize(), x, y, z);
                spawned = solver()->spawnParticle(particle, x, y, z, true);

                if (spawned)
                {
                    ce++;
                }

                c++; //*giggle*
            }

            if (!spawned)
            {
                particle->resetSite();
                delete particle;
                break;
            }

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


uint Boundary::m_maxEventsPrCycle = 1;
uint Boundary::m_coolDown = 1u;
uint Boundary::counter = 0;
