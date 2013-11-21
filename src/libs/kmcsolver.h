#ifndef KMC_SOLVER_H
#define KMC_SOLVER_H

#include <sys/types.h>

#include <armadillo>
using namespace arma;

class Reaction;
class Site;

class KMCSolver
{
public:

    const uint NX;
    const uint NY;
    const uint NZ;

    KMCSolver(uint NX, uint NY, uint NZ);

    void run();

    //REACTION API
    void activateSite(Site *site);
    void deactivateSite(Site *site);

    uint nNeighbours(uint & x, uint & y, uint & z) {
        return neighbours(x, y)(z).n_rows;
    }

    uint nNextNeighbours(uint & x, uint & y, uint & z) {
        return nextNeighbours(x, y)(z).n_rows;
    }




private:

    Site**** sites;

    std::vector<double> allRates;
    std::vector<double> accuAllRates;
    std::vector<uvec> transitions;
    std::vector<uvec> updateRateQueue;
    std::vector<uint> trash;

    ivec delta = {-1, 0, 1};
    ivec delta2 = {-2, -1, 0, 1, 2};

    double t = 0;
    double kTot;

    int counter=0;
    int counter2 = 0;



    uint nTot = 0;
    field<field<umat>> neighbours;
    field<field<umat>> nextNeighbours;
    field<field<umat>> vacantNeighbours;

    void dumpXYZ();

    void getNeighbours(uint i, uint j, uint k);

    void setDiffusionReactions();

    void updateNextNeighbour(uint &x, uint &y, uint &z, const urowvec &newRow, bool activate);


    void updateNeighbourLists(field<field<umat> > &A, field<field<umat> > &B,
                              Site *site, bool activate = false);

    void getRates();

    void updateRates();

    Reaction *getChosenReaction(uint choice);


    void pushToRateQueue(uint &x, uint &y, uint &z);
    void recalcSpecificSite(const uvec &site, uint index);

    bool vectorEqual(const uvec & x, const uvec & y) {
        return ((x(0) == y(0))&&(x(1) == y(1))&&(x(2) == y(2)));
    }
};

#endif // SOLVER_H
