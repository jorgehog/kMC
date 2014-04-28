#include "concentrationwall.h"

#include "../../kmcsolver.h"

#include "../../soluteparticle.h"

#include "../../debugger/debugger.h"

using namespace kMC;

ConcentrationWall::ConcentrationWall(const uint dimension, const uint orientation) :
    Boundary(dimension, orientation, Boundary::ConcentrationWall)
{

}

ConcentrationWall::~ConcentrationWall()
{

}


void ConcentrationWall::update()
{

    KMCDebugger_Assert(m_maxEventsPrCycle, <=, boundarySize(), "Max events pr cycle cannot exceed the number of boundary sites.");


    SoluteParticle *currentParticle;

    uint x, y, z;
    uint c = 0;
    uint ce = 0;

    uint targetN = solver()->targetConcentration()*SoluteParticle::getCurrentSolvantVolume();

    if (SoluteParticle::nSolutionParticles() > targetN)
    {

        while (SoluteParticle::nSolutionParticles() != targetN && c != boundarySize() && ce != m_maxEventsPrCycle)
        {

            getBoundarySite(c, x, y, z);

            currentParticle = solver()->particle(x, y, z);

            if (currentParticle != NULL)
            {
                solver()->despawnParticle(currentParticle);
                ce++;
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

            currentParticle = new SoluteParticle();

            while (c != boundarySize() && !spawned)
            {

                getBoundarySite(c, x, y, z);
                spawned = solver()->spawnParticle(currentParticle, x, y, z, true);

                if (spawned)
                {
                    ce++;
                }

                c++; //*giggle*
            }

            if (!spawned)
            {
                delete currentParticle;
                break;
            }

        }
    }


}




uint ConcentrationWall::m_minDistanceFromSurface;

uint ConcentrationWall::m_maxEventsPrCycle = 1;
