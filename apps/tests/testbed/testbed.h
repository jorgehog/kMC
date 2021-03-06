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

    static void testSaddleOverlapBoxes();

    static void testRateUpdateReach();

    static void testShells();

    static void testTotalParticleStateCounters();

    static void testDistanceTo();

    static void testDiffusionSiteMatrixSetup();

    static void testNeighbors();

    static void testStressedSurface();

    static void testPropertyCalculations();

    static void testRNG();

    static void testParticleMixing();

    static void testBoundarySites();

    static void testConcentrationWall();

    static void testBinarySearchChoise();

    static void testReactionChoise();

    static void testAffectedParticles();

    static void testRateCalculation();

    static void testEnergyAndNeighborSetup();

    static void testUpdateNeigbors();

    static void testInitialReactionSetup();

    static void testSequential();

    static void testKnownCase();

    static void testBoxSizes();

    static void testnNeighborsLimit();

    static void testOptimizedRateVectors();

    static void testReactionVectorUpdate();

    static void testReactionShuffler();

    static void testStateChanges();

    static void testNeighborlist();

    static void testAccuAllRates();

    static void testInitialSiteSetup();


    static void initSimpleSystemParameters(bool clean = true);

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

    static const SnapShot *sequentialCore();

    static void initBoundaryTestParameters();


    static void activateAllSites(const uint sep = 0);

    static void deactivateAllSites();

    static Site * getBoxCenter(const int dx = 0, const int dy = 0, const int dz = 0);

    static SoluteParticle *forceSpawnCenter(const int dx = 0, const int dy = 0, const int dz = 0, const uint particleType = 0);

    static string lastBoundariesName;

    static umat lastBoundaries;

    static void _reactionShufflerCheck(uint nReacs);

    static void forceNewBoxSize(const uvec3 boxSize, bool check = true);

    static void forceNewNNeighborLimit(const uint nNeighborlimit, bool check = true);

    static void forceNewBoundaries(const umat & boundaryMatrix);

    static void forceNewBoundaries(const int boundaryType);

    static double getElectroStaticEnergyContribution(const SoluteParticle *particle);


};
