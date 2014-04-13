#include <kMC>

#include "testbed/testbed.h"

#include <unittest++/UnitTest++.h>
#include <unittest++/TestReporterStdout.h>

#include <iostream>
#include <sys/types.h>


#define focusSuite "Reactions"
//#define focusTest  "ReactionChoise"


#include "defines.h"

SUITE(Misc)
{
    TESTWRAPPER(RNG)

    TESTWRAPPER(BinarySearchChoise)

    TESTWRAPPER(TotalParticleStateCounters)

    TESTWRAPPER(PropertyCalculations)

    TESTWRAPPER(nNeiborsLimit)

    TESTWRAPPER(StateChanges)

    TESTWRAPPER(Neighborlist)

}

SUITE(Reactions)
{

    TESTWRAPPER(RateCalculation)

    TESTWRAPPER(AccuAllRates)

    TESTWRAPPER(ReactionChoise)

    TESTWRAPPER(ReactionVectorUpdate)

    TESTWRAPPER(ReactionShuffler)

}

#define AllBoundaryTests                \
TESTWRAPPER(InitialReactionSetup)       \
                                        \
TESTWRAPPER(BoxSizes)                   \
                                        \
TESTWRAPPER(DistanceTo)                 \
                                        \
TESTWRAPPER(DiffusionSiteMatrixSetup)   \
                                        \
TESTWRAPPER(Neighbors)                  \
                                        \
TESTWRAPPER(UpdateNeigbors)             \
                                        \
TESTWRAPPER(EnergyAndNeighborSetup)     \
                                        \
TESTWRAPPER(OptimizedRateVectors)       \
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

SUITE(ConcWallBoundaries)
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

    mixedBoundaries(2, 0) = Boundary::Edge;
    mixedBoundaries(2, 1) = Boundary::ConcentrationWall;



    int exitSuccess = 0;

    UnitTest::TestReporterStdout reporter;

    UnitTest::TestRunner runner(reporter);

    exitSuccess += RUNSUITE(runner, "Misc");
    exitSuccess += RUNSUITE(runner, "Reactions");

    testBed::initBoundarySuite(Boundary::allBoundariesAs(Boundary::Periodic));
    exitSuccess += RUNSUITE(runner, "PeriodicBoundaries");

    testBed::initBoundarySuite(Boundary::allBoundariesAs(Boundary::Edge));
    exitSuccess += RUNSUITE(runner, "EdgeBoundaries");

    testBed::initBoundarySuite(Boundary::allBoundariesAs(Boundary::ConcentrationWall));
    exitSuccess += RUNSUITE(runner, "ConcWallBoundaries");

    testBed::initBoundarySuite(mixedBoundaries);
    exitSuccess += RUNSUITE(runner, "MixedBoundaries");




    delete testBed::solver;

    return exitSuccess;
}
