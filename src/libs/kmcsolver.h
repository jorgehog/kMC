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

    std::vector<Site*> reactionAffectedSites;

    std::vector<Reaction*> allReactions;

    ivec delta = {-1, 0, 1};

    double t = 0;
    double kTot;

    int counter=0;
    int counter2 = 0;


    void initialize();

    void dumpXYZ();

    void getAllNeighbors();

    void setDiffusionReactions();

    Reaction *getChosenReaction(uint choice);
    uint getReactionChoice(double R);


    void getRateVariables();
    void updateRates();

    bool pushToRateQueue(Site* affectedSite);



};

#endif // SOLVER_H
