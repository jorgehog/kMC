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

    double beta = 0.01;

    uint N = 50;
    uint M = 50;
    uint L = 50;

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
