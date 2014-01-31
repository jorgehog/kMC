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


TEST(NEIGHBOUR_CHECK) {

    testBed test;

    test.testNeighbors();

}

TEST(BINARYSEARCH) {

    testBed test;

    test.testBinarySearchChoise(10000);

}


int main()
{
    return UnitTest::RunAllTests();
}
