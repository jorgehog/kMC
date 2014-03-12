#include <unittest++/UnitTest++.h>

#include <iostream>
#include <sys/types.h>
#include "testbed.h"

#include "../src/libs/debugger/debugger.h"

#include <math.h>

using namespace std;

#define TESTWRAPPER(which) TEST(which) {cout << "Running test " << #which << endl; testBed::test##which();}


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

SUITE(Boundaries)
{
    TEST(Periodic)
    {
        cout << "Running test Periodic" << endl;
        testBed::runAllBoundaryTests(zeros<umat>(3, 2) + Boundary::Periodic);
    }

    TEST(Edge)
    {
        cout << "Running test Edge" << endl;
        testBed::runAllBoundaryTests(zeros<umat>(3, 2) + Boundary::Edge);
    }

    TEST(Surface)
    {
        cout << "Running test Surface" << endl;
        testBed::runAllBoundaryTests(zeros<umat>(3, 2) + Boundary::Surface);
    }

    TEST(Mixed)
    {
        cout << "Running test Mixed" << endl;
        umat mixedBoundaries(3, 2);

        mixedBoundaries(0, 0) = Boundary::Periodic;
        mixedBoundaries(0, 1) = Boundary::Periodic;

        mixedBoundaries(1, 0) = Boundary::Edge;
        mixedBoundaries(1, 1) = Boundary::Edge;

        mixedBoundaries(2, 0) = Boundary::Surface;
        mixedBoundaries(2, 1) = Boundary::ConcentrationWall;

        testBed::runAllBoundaryTests(mixedBoundaries);

    }

}

SUITE(General)
{
    TESTWRAPPER(Sequential)

    TESTWRAPPER(KnownCase)
}

int main()
{
    KMCDebugger_SetFilename("testTrace");
    return UnitTest::RunAllTests();
}
