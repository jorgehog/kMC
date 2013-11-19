#ifndef REACTION_H
#define REACTION_H

#include <armadillo>
using namespace arma;

class Site;
class KMCSolver;

class Reaction
{
public:
    Reaction(Site* site);

    virtual double calcRate() = 0;

    void setKMCSolverObject(KMCSolver * solver) {
        this->solver = solver;
    }

private:

    KMCSolver* solver;

};

#endif // REACTION_H
