#include <unittest++/UnitTest++.h>

#include <iostream>
#include <sys/types.h>
#include "testbed.h"

#include "../../src/libs/debugger/debugger.h"

#include <math.h>

using namespace std;

#define UBERTEST(which)                                 \
TEST(which)                                             \
{                                                       \
    cout << "Testing: " << #which << endl;              \
    testBed test;                                       \
    test.test##which();                                 \
    cout << "------------------------------\n" << endl; \
}

UBERTEST(RNG)
UBERTEST(EnergyAndNeighborSetup)
UBERTEST(DiffusionSiteMatrixSetup)
UBERTEST(Neighbors)
UBERTEST(DistanceTo)
UBERTEST(BinarySearchChoise)
UBERTEST(UpdateNeigbors)
UBERTEST(RateCalculation)
UBERTEST(ReactionChoise)
UBERTEST(HasCrystalNeighbor)
UBERTEST(InitializationOfCrystal)
UBERTEST(InitialReactionSetup)
UBERTEST(Sequential)
UBERTEST(KnownCase)

int main()
{
    KMCDebugger_SetFilename("testTrace");
    return UnitTest::RunAllTests();
}
