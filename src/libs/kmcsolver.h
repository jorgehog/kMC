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

    void addReaction(Reaction* reaction) {
        allReactions.push_back(reaction);
    }

    friend class testBed;


private:

    Site**** sites;

    std::vector<double> accuAllRates;

    std::vector<Site*> reactionAffectedSites;

    std::vector<Reaction*> allReactions;

    ivec delta = {-1, 0, 1};

    double t = 0;
    double kTot;

    int counter=0;
    int counter2 = 0;

    uint nTot = 0;
    field<field<umat>> neighbours;
    field<field<umat>> nextNeighbours;
    field<field<umat>> vacantNeighbours;

    void initialize();

    void dumpXYZ();

    void getAllNeighbours();
    void getNeighbours(uint i, uint j, uint k);

    void updateNextNeighbour(uint &x, uint &y, uint &z, const urowvec &newRow, bool activate);

    void updateNeighbourLists(field<field<umat> > &A, field<field<umat> > &B,
                              Site *site, bool activate = false);


    void setDiffusionReactions();

    Reaction *getChosenReaction(uint choice);
    uint getReactionChoice(double R);


    void getRateVariables();
    void updateRates();

    bool pushToRateQueue(Site* affectedSite);

};

#endif // SOLVER_H
