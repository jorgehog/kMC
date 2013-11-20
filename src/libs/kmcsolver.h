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

    KMCSolver(uint NX, uint NY, uint NZ);

    void addReaction(Reaction* reaction, uint & x, uint &y, uint &z);

    Site**** getSites() const {
        return sites;
    }


    void test();

    void run();

    //REACTION API
    void activateSite(uint i, uint j, uint k);
    void deactivateSite(uint i, uint j, uint k);

private:

    uint NX;
    uint NY;
    uint NZ;

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

    double Enn = 2;
    double Ennn = 1.7;
    double EspN = 1;
    double EspNN = 0;
    double temperature = 10;
    double mu = 1;

    uint nTot = 0;
    field<field<umat>> neighbours;
    field<field<umat>> nextNeighbours;
    field<field<umat>> vacantNeighbours;

    void dumpXYZ();

    void getNeighbours(uint i, uint j, uint k);
    void updateNextNeighbour(uint &x, uint &y, uint &z, const urowvec &newRow, bool activate);

    void reactionDiffusion(uint i, uint j, uint k);

    void updateNeighbourLists(field<field<umat> > &A, field<field<umat> > &B,
                              uint i, uint j, uint k, bool activate = false);

    void getAllTransitionsAndRates();
    void getAllTransitionsAndRatesForSite(uint i, uint j, uint k);
    void updateTransitionsAndRates();
    void pushToRateQueue(uint &x, uint &y, uint &z);
    void recalcSpecificSite(const uvec &site, uint index);

    bool vectorEqual(const uvec & x, const uvec & y) {
        return ((x(0) == y(0))&&(x(1) == y(1))&&(x(2) == y(2)));
    }
};

#endif // SOLVER_H
