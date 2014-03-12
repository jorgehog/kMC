#pragma once


#include <kMC>

#include <sys/types.h>

#define MIN(x, y) x < y ? x : y

#include <libconfig.h++>

using namespace libconfig;
using namespace kMC;

class testBed
{
public:

    static void makeSolver();

    static void testDistanceTo();

    static void testDeactivateSurface();

    static void testDiffusionSiteMatrixSetup();

    static void testNeighbors();

    static void testRNG();

    static void testBinarySearchChoise();

    static void testReactionChoise();

    static void testRateCalculation();

    static void testEnergyAndNeighborSetup();

    static void testUpdateNeigbors();

    static void testHasCrystalNeighbor();

    static void testInitializationOfCrystal();

    static void testInitialReactionSetup();

    static void testSequential();

    static void testKnownCase();

    static void testBoxSizes();

    static void testnNeiborsLimit();

    static void testnNeighborsToCrystallize();

    static void testBoundaries();

    static void testDiffusionSeparation();

    static void runAllBoundaryTests();

    static KMCSolver* solver;


    static const inline uint & NX()
    {
        return solver->NX();
    }

    static const inline  uint & NY()
    {
        return solver->NY();
    }

    static const inline uint & NZ()
    {
        return solver->NZ();
    }

};
