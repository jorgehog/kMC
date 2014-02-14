#include <unittest++/UnitTest++.h>

#include <iostream>
#include <sys/types.h>
#include "testbed.h"

#include <math.h>


using namespace std;


TEST(RNG_CHECK) {

    testBed test;

    test.testRNG();

}

TEST(NEIGHBORS_SETUP) {
    testBed test;
    test.testEnergyAndNeighborSetup();
}

TEST(NEIGHBOUR_CHECK) {

    testBed test;

    test.testNeighbors();

}

TEST(DISTANCE_TO) {
    testBed test;

    test.testDistanceTo();

}

TEST(BINARYSEARCH) {

    testBed test;

    test.testBinarySearchChoise(10000);

}

TEST(UPDATE_NEIGHBORS) {
    testBed test;
    test.testUpdateNeigbors();
}

TEST(RATECALC) {
    testBed test;
    test.testRateCalculation();
}

TEST(REACTIONCHOISE) {

    testBed test;

    test.testReactionChoise(1);
}

TEST(CRYSTALSTUFF1) {
    testBed test;
    test.testHasCrystalNeighbor();
}

TEST(INITCRYSTAL) {
    testBed test;
    test.testInitializationOfCrystal();
}

TEST(INITREACTIONS)
{
    testBed test;
    test.testInitialReactionSetup();
}


TEST(KnownCase)
{
    testBed test;
    test.testKnownCase();
}

int main()
{
    return UnitTest::RunAllTests();
}
