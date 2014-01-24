#ifndef KMC_SOLVER_H
#define KMC_SOLVER_H

#include <sys/types.h>
#include <armadillo>

#include "site.h"

using namespace arma;

class Reaction;

class KMCSolver
{
public:

    const uint NX;
    const uint NY;
    const uint NZ;


    KMCSolver(uint NX, uint NY, uint NZ);


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

    friend class testBed;


private:

    const uint nNeighborLimit = 2;


    Site**** sites;


    std::vector<double> accuAllRates;

    std::vector<Reaction*> allReactions;


    double kTot;
    double totalTime = 0;

    uint nCycles = 100000;
    uint cycle   = 0;

    uint cyclesPerOutput = 100;
    uint outputCounter   = 0;


    void initialize();

    void getAllNeighbors();

    void setDiffusionReactions();

    void getRateVariables();

    uint getReactionChoice(double R);

    Reaction *getChosenReaction(uint choice);

    void dumpXYZ();

    void dumpOutput();


};

#endif // SOLVER_H
