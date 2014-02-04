#ifndef TESTBED_H
#define TESTBED_H

#include <sys/types.h>

class KMCSolver;

class testBed
{
public:
    testBed();

    ~testBed();

    void testNeighbors();

    void testRNG();

    void testBinarySearchChoise(uint LIM);


    void testReactionChoise(uint LIM);

    void testRateCalculation();

    void testEnergyAndNeighborSetup();

    void testUpdateNeigbors();



    uint failCount;
    uint winCount;
    uint nTrials;

    void reset() {

        failCount = 0;
        winCount = 0;
        nTrials = 0;
    }

    uint NX;
    uint NY;
    uint NZ;

    KMCSolver* solver;
};

#endif // TESTBED_H
