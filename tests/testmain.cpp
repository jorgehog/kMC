#include <unittest++/UnitTest++.h>

#include <iostream>
#include <sys/types.h>
#include "testbed.h"

#include "../src/libs/debugger/debugger.h"

#include <math.h>

using namespace std;

#define TESTWRAPPER(which)                                 \
TEST(which)                                             \
{                                                       \
    cout << "Testing: " << #which << endl;              \
    testBed test;                                       \
    test.test##which();                                 \
    cout << "------------------------------\n" << endl; \
}

TESTWRAPPER(DeactivateSurface)
TESTWRAPPER(RNG)
TESTWRAPPER(EnergyAndNeighborSetup)
TESTWRAPPER(DiffusionSiteMatrixSetup)
TESTWRAPPER(Neighbors)
TESTWRAPPER(DistanceTo)
TESTWRAPPER(BinarySearchChoise)
TESTWRAPPER(UpdateNeigbors)
TESTWRAPPER(RateCalculation)
TESTWRAPPER(ReactionChoise)
TESTWRAPPER(HasCrystalNeighbor)
TESTWRAPPER(InitializationOfCrystal)
TESTWRAPPER(InitialReactionSetup)
TESTWRAPPER(Sequential)
TESTWRAPPER(BoxSizes)
TESTWRAPPER(nNeiborsLimit)
TESTWRAPPER(nNeighborsToCrystallize)
TESTWRAPPER(Boundaries)
TESTWRAPPER(DiffusionSeparation)
TESTWRAPPER(KnownCase)

int main()
{
    KMCDebugger_SetFilename("testTrace");
    return UnitTest::RunAllTests();
}
