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


TEST(NEIGHBOURUPDATE_CHECK) {

    testBed test;

    test.testNeighboursUpdating();



}

TEST(NEIGHBOUR_CHECK) {

    testBed test;

    test.testNeighbours();

}


int main()
{
    return UnitTest::RunAllTests();
}
