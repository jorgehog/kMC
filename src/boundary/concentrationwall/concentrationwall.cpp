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

    KMCDebugger_Assert(m_maxEventsPrCycle, <=, boundarySites().size(), "Max events pr cycle cannot exceed the number of boundary sites.");


    Site * currentSite;

    uint c = 0;
    uint ce = 0;


    //    crystalBoxTopology = Site::getCurrentCrystalBoxTopology();

    //    bool resize;

    //    switch (orientation()) {
    //    case Near:
    //        resize = crystalBoxTopology(dimension(), Near) < minDistanceFromSurface;

    //        break;
    //    case Far:
    //        resize = crystalBoxTopology(dimension(), Far) > span() - minDistanceFromSurface;

    //        break;
    //    }


    //    if (resize)
    //    {

    //        uvec3 N = solver()->NVec();

    //        N(dimension()) += systemSizeIncrementSize; //Size size size...

    //        solver()->setBoxSize(N, true, true);

    //    }


    std::random_shuffle(boundarySites().begin(), boundarySites().end(), [] (uint n) {return KMC_RNG_UNIFORM()*n;});

    if (SoluteParticle::getCurrentConcentration() > solver()->targetConcentration())
    {

        while (SoluteParticle::getCurrentConcentration() > solver()->targetConcentration() && c != boundarySites().size() && ce != m_maxEventsPrCycle)
        {
            currentSite = boundarySites().at(c);

            if (currentSite->isActive())
            {
                solver()->despawnParticle(currentSite);
                ce++;
            }

            c++; //*giggle*
        }

    }

    else
    {

        bool spawned;

        while (ce != m_maxEventsPrCycle && SoluteParticle::getCurrentConcentration() < solver()->targetConcentration())
        {

            spawned = false;

            SoluteParticle *particle = new SoluteParticle();

            while (c != boundarySites().size() && !spawned)
            {
                currentSite = boundarySites().at(c);

                spawned = solver()->spawnParticle(particle, currentSite, true);

                if (spawned)
                {
                    ce++;
                }

                c++; //*giggle*
            }

            if (!spawned)
            {
                delete particle;
                break;
            }

        }
    }
}

void ConcentrationWall::initialize()
{
    setupBoundarySites();
}



uint ConcentrationWall::m_minDistanceFromSurface;

uint ConcentrationWall::m_maxEventsPrCycle = 3;
