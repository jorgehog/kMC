#include <unittest++/UnitTest++.h>

#include <iostream>
#include <sys/types.h>
#include "testbed.h"

#include "../src/libs/debugger/debugger.h"




//#define focusSuite "Misc"
//#define focusTest  "PropertyCalculations"

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

SUITE(PeriodicBoundaries)
{
    TESTWRAPPER(RunAllBoundaryTests, zeros<umat>(3, 2) + Boundary::Periodic)
}

SUITE(EdgeBoundaries)
{
    TESTWRAPPER(RunAllBoundaryTests, zeros<umat>(3, 2) + Boundary::Edge)
}

SUITE(ConcentrationBoundaries)
{
    TESTWRAPPER(RunAllBoundaryTests, zeros<umat>(3, 2) + Boundary::ConcentrationWall)
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


int main()
{
    using namespace SuiteMixedBoundaries;

    KMCDebugger_SetFilename("testTrace");


    mixedBoundaries(0, 0) = Boundary::Periodic;
    mixedBoundaries(0, 1) = Boundary::Periodic;

    mixedBoundaries(1, 0) = Boundary::Edge;
    mixedBoundaries(1, 1) = Boundary::Edge;

    mixedBoundaries(2, 0) = Boundary::Surface;
    mixedBoundaries(2, 1) = Boundary::ConcentrationWall;


    testBed::makeSolver();

    int exitSuccess = UnitTest::RunAllTests();

    delete testBed::solver;


    return exitSuccess;
}
