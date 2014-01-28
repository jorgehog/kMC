#include <unittest++/UnitTest++.h>

#include <iostream>
#include <sys/types.h>
#include "testbed.h"

#include <math.h>

#ifndef UNITTESTCPP_H
#define UNITTESTCPP_H
#define TEST(name) \
    void lol#name
#endif


using namespace std;

TEST(RNG_CHECK) {

    testBed test;

    test.testRNG();

}


TEST(NEIGHBOUR_CHECK) {

    testBed test;

    test.testNeighbors();

}


int main()
{
    return UnitTest::RunAllTests();
}
