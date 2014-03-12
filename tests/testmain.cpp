#include <unittest++/UnitTest++.h>

#include <iostream>
#include <sys/types.h>
#include "testbed.h"

#include "../src/libs/debugger/debugger.h"

#include <math.h>

using namespace std;

SUITE(Misc)
{

    TEST(RNG)
    {
        testBed::testRNG();
    }

    TEST(BinarySearchChoise)
    {
        testBed::testBinarySearchChoise();
    }

}

SUITE(Reactions)
{
    TEST(RateCalculation)
    {
        testBed::testRateCalculation();
    }

    TEST(ReactionChoise)
    {
        testBed::testReactionChoise();
    }

    TEST(InitialReactionSetup)
    {
        testBed::testInitialReactionSetup();
    }
}

SUITE(StateChanges)
{
    TEST(HasCrystalNeighbor)
    {
        testBed::testHasCrystalNeighbor();
    }

    TEST(DeactivateSurface)
    {
        testBed::testDeactivateSurface();
    }
}

SUITE(Parameters)
{
    TEST(BoxSizes)
    {
        testBed::testBoxSizes();
    }

    TEST(nNeiborsLimit)
    {
        testBed::testnNeiborsLimit();
    }

    TEST(nNeighborsToCrystallize)
    {
        testBed::testnNeighborsToCrystallize();
    }

    TEST(DiffusionSeparation)
    {
        testBed::testDiffusionSeparation();
    }
}

SUITE(Boundaries)
{
    TEST(Periodic)
    {
        Site::setBoundaries(zeros<umat>(3, 2) + Boundary::Periodic);
        testBed::runAllBoundaryTests();
    }

    TEST(Edge)
    {
            Site::setBoundaries(zeros<umat>(3, 2) + Boundary::Edge);
        testBed::runAllBoundaryTests();
    }

    TEST(Surface)
    {
            Site::setBoundaries(zeros<umat>(3, 2) + Boundary::Surface);
        testBed::runAllBoundaryTests();
    }

    TEST(Mixed)
    {
        umat mixedBoundaries(3, 2);

        mixedBoundaries(0, 0) = Boundary::Periodic;
        mixedBoundaries(0, 1) = Boundary::Periodic;

        mixedBoundaries(1, 0) = Boundary::Edge;
        mixedBoundaries(1, 1) = Boundary::Edge;

        mixedBoundaries(2, 0) = Boundary::Surface;
        mixedBoundaries(2, 1) = Boundary::ConcentrationWall;

        Site::setBoundaries(mixedBoundaries);

        testBed::runAllBoundaryTests();

    }

}

SUITE(General)
{
    TEST(Sequential)
    {
        testBed::testSequential();
    }

    TEST(KnownCase)
    {
        testBed::testKnownCase();
    }
}

int main()
{
    KMCDebugger_SetFilename("testTrace");
    return UnitTest::RunAllTests();
}
