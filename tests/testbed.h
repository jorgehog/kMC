#ifndef TESTBED_H
#define TESTBED_H

#include <sys/types.h>

class testBed
{
public:
    testBed();

    void testNeighbors();

    void testRNG();

    uint failCount;
    uint winCount;
    uint nTrials;

    void reset() {

        failCount = 0;
        winCount = 0;
        nTrials = 0;
    }
};

#endif // TESTBED_H
