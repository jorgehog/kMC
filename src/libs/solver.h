#ifndef KMC_SOLVER_H
#define KMC_SOLVER_H

#include <sys/types.h>

#include <armadillo>
using namespace arma;

typedef bool*** BoolCube;

class Solver
{
public:

    Solver();

    double alpha = 1;
    double beta = 0.1;

    uint N = 30;
    uint M = 30;
    uint L = 30;

    ivec delta = {-1, 0, 1};

    double t = 0;
    int counter=0;
    int counter2 = 0;

    uint nTot = 0;

    BoolCube siteCube;
    field<field<umat>> neighbours;
    field<field<umat>> vacantNeighbours;

    void test();
    void dump();

    double getRate(double E);

    void getNeighbours(uint i, uint j, uint k);

    void reactionDiffusion(uint i, uint j, uint k);

    void reactionCreation(uint i, uint j, uint k);
    void reactionDeletion(uint i, uint j, uint k);

    void updateNeighbourLists(field<field<umat> > & A, field<field<umat> > & B, uint i, uint j, uint k);
};

#endif // SOLVER_H
