#include <unittest++/UnitTest++.h>

#include <iostream>
#include <sys/types.h>
#include "testbed.h"

#include "../src/libs/debugger/debugger.h"

#include <math.h>

using namespace std;


//#define focustest "nNeighborsToCrystallize"


#ifdef focustest
#define shouldrun(which) (strcmp(#which, focustest) == 0)
#else
#define shouldrun(which) true
#endif

#define TESTWRAPPER(which)                              \
TEST(which)                                             \
{                                                       \
    if (!shouldrun(which)) return;                      \
    \
    cout << "Testing: " << #which << endl;              \
    testBed test;                                       \
    test.test##which();                                 \
    cout << "------------------------------\n" << endl; \
}

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
TESTWRAPPER(DeactivateSurface)
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
