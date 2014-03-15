#pragma once


#include <kMC>

#include <sys/types.h>

#define MIN(x, y) x < y ? x : y

#include <libconfig.h++>

using namespace libconfig;
using namespace kMC;

class SnapShot;

class testBed
{
public:

    static void makeSolver();

    static void testTotalParticleStateCounters();

    static void testDistanceTo();

    static void testDeactivateSurface();

    static void testDiffusionSiteMatrixSetup();

    static void testNeighbors();

    static void testPropertyCalculations();

    static void testRNG();

    static void testBinarySearchChoise();

    static void testReactionChoise();

    static void testRateCalculation();

    static void testEnergyAndNeighborSetup();

    static void testUpdateNeigbors();

    static void testHasCrystalNeighbor();

    static void testInitializationOfCrystal();

    static void testInitialReactionSetup();

    static void testSequential(const umat &boundaries);

    static void testKnownCase(const umat &boundaries, const string name);

    static void testBoxSizes();

    static void testnNeiborsLimit();

    static void testnNeighborsToCrystallize();

    static void testDiffusionSeparation();

    static void testRunAllBoundaryTests(const umat &boundaries);

    static KMCSolver* solver;

    static wall_clock timer;

    static seed_type baseSeed;


private:

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

    static const SnapShot *testSequentialCore();

    static void initBoundaryTestParameters(const umat &boundaries);

    static void initSimpleSystemParameters();

    static void activateAllSites();

    static void deactivateAllSites();

    static Site * getBoxCenter();


};
