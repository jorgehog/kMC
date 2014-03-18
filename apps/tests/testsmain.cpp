#include <kMC>

#include "testbed.h"

#include <unittest++/UnitTest++.h>
#include <unittest++/TestReporterStdout.h>

#include <iostream>
#include <sys/types.h>


//#define focusSuite "SurfaceBoundaries"
//#define focusTest  "KnownCase"


#include "defines.h"

SUITE(Misc)
{
    TESTWRAPPER(RNG)

    TESTWRAPPER(BinarySearchChoise)

    TESTWRAPPER(TotalParticleStateCounters)

    TESTWRAPPER(PropertyCalculations)
}

SUITE(Reactions)
{

    TESTWRAPPER(RateCalculation)

    TESTWRAPPER(ReactionChoise)

}

SUITE(StateChanges)
{
    TESTWRAPPER(HasCrystalNeighbor)

    TESTWRAPPER(DeactivateSurface)
}

SUITE(Parameters)
{

    TESTWRAPPER(nNeiborsLimit)

    TESTWRAPPER(nNeighborsToCrystallize)

    TESTWRAPPER(DiffusionSeparation)

}

#define AllBoundaryTests                \
TESTWRAPPER(InitialReactionSetup)       \
                                        \
TESTWRAPPER(BoxSizes)                   \
                                        \
TESTWRAPPER(DistanceTo)                 \
                                        \
TESTWRAPPER(InitializationOfCrystal)    \
                                        \
TESTWRAPPER(DiffusionSiteMatrixSetup)   \
                                        \
TESTWRAPPER(Neighbors)                  \
                                        \
TESTWRAPPER(UpdateNeigbors)             \
                                        \
TESTWRAPPER(EnergyAndNeighborSetup)     \
                                        \
TESTWRAPPER(Sequential)                 \
                                        \
TESTWRAPPER(KnownCase)


SUITE(PeriodicBoundaries)
{
    AllBoundaryTests
}

SUITE(EdgeBoundaries)
{
    AllBoundaryTests
}

SUITE(ConcentrationBoundaries)
{
    AllBoundaryTests
}

SUITE(SurfaceBoundaries)
{
    AllBoundaryTests
}

SUITE(MixedBoundaries)
{
    AllBoundaryTests
}

#define AllBoundariesAs(type) (zeros<umat>(3, 2) + type)

int main()
{
    using namespace SuiteMixedBoundaries;

    KMCDebugger_SetEnabledTo(false);

    testBed::makeSolver();

    umat mixedBoundaries(3, 2);

    mixedBoundaries(0, 0) = Boundary::Periodic;
    mixedBoundaries(0, 1) = Boundary::Periodic;

    mixedBoundaries(1, 0) = Boundary::Edge;
    mixedBoundaries(1, 1) = Boundary::Edge;

    mixedBoundaries(2, 0) = Boundary::Surface;
    mixedBoundaries(2, 1) = Boundary::ConcentrationWall;



    int exitSuccess = 0;

    UnitTest::TestReporterStdout reporter;

    UnitTest::TestRunner runner(reporter);

    exitSuccess += RUNSUITE(runner, "Misc");
    exitSuccess += RUNSUITE(runner, "Reactions");
    exitSuccess += RUNSUITE(runner, "StateChanges");
    exitSuccess += RUNSUITE(runner, "Parameters");

    testBed::initBoundarySuite(AllBoundariesAs(Boundary::Periodic));
    exitSuccess += RUNSUITE(runner, "PeriodicBoundaries");

    testBed::initBoundarySuite(AllBoundariesAs(Boundary::Edge));
    exitSuccess += RUNSUITE(runner, "EdgeBoundaries");

    testBed::initBoundarySuite(AllBoundariesAs(Boundary::Surface));
    exitSuccess += RUNSUITE(runner, "SurfaceBoundaries");

    testBed::initBoundarySuite(AllBoundariesAs(Boundary::ConcentrationWall));
    exitSuccess += RUNSUITE(runner, "ConcentrationBoundaries");

    testBed::initBoundarySuite(mixedBoundaries);
    exitSuccess += RUNSUITE(runner, "MixedBoundaries");




    delete testBed::solver;

    return exitSuccess;
}
