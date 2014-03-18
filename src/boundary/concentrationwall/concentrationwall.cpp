#include "concentrationwall.h"

#include "../../kmcsolver.h"


using namespace kMC;

ConcentrationWall::ConcentrationWall(const uint dimension, const uint orientation) :
    Boundary(dimension, orientation, Boundary::ConcentrationWall)
{

}

ConcentrationWall::~ConcentrationWall()
{

}

void ConcentrationWall::loadConfig(const Setting &setting)
{
//    try
//    {

//        minDistanceFromSurface = setting["ds"];
//    }
//    catch (const SettingNotFoundException & exc)
//    {

//        minDistanceFromSurface = span()/4;
//        return;

//    }
}

void ConcentrationWall::update()
{

    Site * currentSite;

    uint c = 0;


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


    while (Site::getCurrentSolutionDensity() > solver()->targetSaturation() && c != boundarySites().size())
    {
        currentSite = boundarySites().at(c);

        if (currentSite->isActive())
        {
            currentSite->deactivate();
        }

        c++; //*giggle*
    }

    while (Site::getCurrentSolutionDensity() < solver()->targetSaturation() && c != boundarySites().size())
    {
        currentSite = boundarySites().at(c);

        if (currentSite->isLegalToSpawn())
        {
            currentSite->activate();
        }

        c++; //*giggle*
    }

}

void ConcentrationWall::initialize()
{
    setupBoundarySites();
}