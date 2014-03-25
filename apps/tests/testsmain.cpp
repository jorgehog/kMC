#include <kMC>

#include "testbed/testbed.h"

#include <unittest++/UnitTest++.h>
#include <unittest++/TestReporterStdout.h>

#include <iostream>
#include <sys/types.h>


#define focusSuite "PeriodicBoundaries"
#define focusTest  "KnownCase"


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

SUITE(ConcentrationWallBoundaries)
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

    testBed::initBoundarySuite(Boundary::allBoundariesAs(Boundary::Periodic));
    exitSuccess += RUNSUITE(runner, "PeriodicBoundaries");

    testBed::initBoundarySuite(Boundary::allBoundariesAs(Boundary::Edge));
    exitSuccess += RUNSUITE(runner, "EdgeBoundaries");

    testBed::initBoundarySuite(Boundary::allBoundariesAs(Boundary::Surface));
    exitSuccess += RUNSUITE(runner, "SurfaceBoundaries");

    testBed::initBoundarySuite(Boundary::allBoundariesAs(Boundary::ConcentrationWall));
    exitSuccess += RUNSUITE(runner, "ConcentrationWallBoundaries");

    testBed::initBoundarySuite(mixedBoundaries);
    exitSuccess += RUNSUITE(runner, "MixedBoundaries");




    delete testBed::solver;

    return exitSuccess;
}
