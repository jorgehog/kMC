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

    static void testSequential();

    static void testKnownCase();

    static void testBoxSizes();

    static void testnNeiborsLimit();

    static void testnNeighborsToCrystallize();

    static void testDiffusionSeparation();

    static void testOptimizedRateVectors();

    static void testReactionVectorUpdate();

    static void testReactionShuffler();




    static void initBoundarySuite(const umat &boundaries);

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

    static void mainloop_meat();

    static void fill_rate_stuff(vector<double> &accuAllRates, vector<Reaction *> &allPossibleReactions, double &kTot);

    static const SnapShot *testSequentialCore();

    static void initBoundaryTestParameters();

    static void initSimpleSystemParameters();

    static void activateAllSites();

    static void deactivateAllSites();

    static Site * getBoxCenter(const int dx = 0, const int dy = 0, const int dz = 0);

    static string lastBoundariesName;

    static umat lastBoundaries;

    static void _reactionShufflerCheck(uint nReacs);

};
