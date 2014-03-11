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
    testBed();

    KMCSolver * makeSolver();

    ~testBed();

    void testDistanceTo();

    void testDeactivateSurface();

    void testDiffusionSiteMatrixSetup();

    void testNeighbors();

    void testRNG();

    void testBinarySearchChoise();

    void testReactionChoise();

    void testRateCalculation();

    void testEnergyAndNeighborSetup();

    void testUpdateNeigbors();

    void testHasCrystalNeighbor();

    void testInitializationOfCrystal();

    void testInitialReactionSetup();

    void testSequential();

    void testKnownCase();

    void testBoxSizes();

    void testnNeiborsLimit();

    void testnNeighborsToCrystallize();

    void testPeriodicBoundaries();

    void testOpenBoundaries();

    void testWallBoundaries();

    void testConcentrationBoundaries();

    void testMixedboundaries();

    void testDiffusionSeparation();





    uint failCount;
    uint winCount;
    uint nTrials;

    void reset() {

        failCount = 0;
        winCount = 0;
        nTrials = 0;
    }

    const uint &NX () const
    {
        return solver->NX();
    }

    const uint &NY () const
    {
        return solver->NY();
    }

    const uint &NZ () const
    {
        return solver->NZ();
    }

    KMCSolver* solver;
    Setting* root;
};
