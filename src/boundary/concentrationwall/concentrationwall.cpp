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
    counter++;
    if (counter%m_coolDown != 0)
    {
        return;
    }

    KMCDebugger_Assert(m_maxEventsPrCycle, <=, boundarySize(), "Max events pr cycle cannot exceed the number of boundary sites.");


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

            getBoundarySite(c, x, y, z);

            currentSite = solver()->getSite(x, y, z);

            if (currentSite->isActive())
            {
                solver()->despawnParticle(currentSite->associatedParticle());
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

            SoluteParticle *particle = new SoluteParticle();

            while (c != boundarySize() && !spawned)
            {

                getBoundarySite(c, x, y, z);
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

uint ConcentrationWall::m_maxEventsPrCycle = 1;
uint ConcentrationWall::counter = 0;
