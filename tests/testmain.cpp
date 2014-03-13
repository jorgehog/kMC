#include <unittest++/UnitTest++.h>

#include <iostream>
#include <sys/types.h>
#include "testbed.h"

#include "../src/libs/debugger/debugger.h"

#include <math.h>

using namespace std;

#define TESTCORE(which, ...) \
    cout << "Running test " << #which << endl; \
    \
    testBed::makeSolver(); \
    \
    testBed::timer.tic(); \
    testBed::test##which(__VA_ARGS__); \
    cout << "Done (" << testBed::timer.toc() << " s)" << endl; \
    \
    delete testBed::solver; \
    \


//Defined in one line to made unittest++ file line match.
#define TESTWRAPPER(which, ...) TEST(which) {TESTCORE(which, ##__VA_ARGS__)}


SUITE(Misc)
{
    TESTWRAPPER(RNG)

    TESTWRAPPER(BinarySearchChoise)
}

SUITE(Reactions)
{
    TESTWRAPPER(RateCalculation)

    TESTWRAPPER(ReactionChoise)

    TESTWRAPPER(InitialReactionSetup)
}

SUITE(StateChanges)
{
    TESTWRAPPER(HasCrystalNeighbor)

    TESTWRAPPER(DeactivateSurface)
}

SUITE(Parameters)
{
    TESTWRAPPER(BoxSizes)

    TESTWRAPPER(nNeiborsLimit)

    TESTWRAPPER(nNeighborsToCrystallize)

    TESTWRAPPER(DiffusionSeparation)

}

SUITE(PeriodicBoundaries)
{
    TESTWRAPPER(RunAllBoundaryTests, zeros<umat>(3, 2) + Boundary::Periodic)
}

SUITE(EdgeBoundaries)
{
    TESTWRAPPER(RunAllBoundaryTests, zeros<umat>(3, 2) + Boundary::Edge)
}
SUITE(SurfaceBoundaries)
{
    TESTWRAPPER(RunAllBoundaryTests, zeros<umat>(3, 2) + Boundary::Surface)
}

SUITE(MixedBoundaries)
{
    umat mixedBoundaries(3, 2);
    TESTWRAPPER(RunAllBoundaryTests, mixedBoundaries)
}


SUITE(General)
{
    TESTWRAPPER(Sequential)

    TESTWRAPPER(KnownCase)
}

int main()
{
    KMCDebugger_SetFilename("testTrace");

    using namespace SuiteMixedBoundaries;

    mixedBoundaries(0, 0) = Boundary::Periodic;
    mixedBoundaries(0, 1) = Boundary::Periodic;

    mixedBoundaries(1, 0) = Boundary::Edge;
    mixedBoundaries(1, 1) = Boundary::Edge;

    mixedBoundaries(2, 0) = Boundary::Surface;
    mixedBoundaries(2, 1) = Boundary::ConcentrationWall;

    return UnitTest::RunAllTests();
}
