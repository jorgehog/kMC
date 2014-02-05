#ifndef KMC_SOLVER_H
#define KMC_SOLVER_H

#include <sys/types.h>
#include <armadillo>

#include <libconfig_utils/libconfig_utils.h>

#include "site.h"

using namespace arma;

class Reaction;

class KMCSolver
{
public:


    KMCSolver(const Setting & root);

    ~KMCSolver();

    void run();


    uint nNeighbors(uint & x, uint & y, uint & z) {
        return sites[x][y][z]->nNeighbors(0);
    }

    uint nNextNeighbors(uint & x, uint & y, uint & z) {
        return sites[x][y][z]->nNeighbors(1);
    }


    void addReaction(Reaction* reaction) {
        allReactions.push_back(reaction);
    }


    Site**** getSites() {
        return sites;
    }

    uint getNX () {
        return NX;
    }

    uint getNY () {
        return NY;
    }

    uint getNZ () {
        return NZ;
    }

    friend class testBed;


private:

    uint NX;
    uint NY;
    uint NZ;

    Site**** sites;


    std::vector<double> accuAllRates;

    std::vector<Reaction*> allReactions;


    double kTot;
    double totalTime = 0;

    uint nCycles;
    uint cycle   = 0;

    uint cyclesPerOutput;
    uint outputCounter   = 0;

    double saturation;
    double RelativeSeedSize;

    void initializeCrystal();

    void spawnCrystalSeed();

    void getAllNeighbors();

    void setDiffusionReactions();

    void getRateVariables();

    uint getReactionChoice(double R);



    void dumpXYZ();

    void dumpOutput();


};

#endif // SOLVER_H
