#pragma once


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


    uint nNeighbors(uint & x, uint & y, uint & z)
    {
        return sites[x][y][z]->nNeighbors(0);
    }

    uint nNextNeighbors(uint & x, uint & y, uint & z)
    {
        return sites[x][y][z]->nNeighbors(1);
    }

    Site**** getSites()
    {
        return sites;
    }

    const uint &getNX ()
    {
        return NX;
    }

    const uint &getNY ()
    {
        return NY;
    }

    const uint &getNZ ()
    {
        return NZ;
    }

    friend class testBed;


private:

    double saturation;

    double RelativeSeedSize;


    Site**** sites;

    uint NX;
    uint NY;
    uint NZ;


    double kTot;
    std::vector<double> accuAllRates;


    std::vector<Reaction*> allReactions;


    double totalTime;

    uint nCycles;
    uint cycle;

    uint cyclesPerOutput;
    uint outputCounter;


    void initializeCrystal();

    void initializeDiffusionReactions();

    void initializeSites();


    void getRateVariables();

    uint getReactionChoice(double R);



    void dumpXYZ();

    void dumpOutput();

    static uint ptrCount;

};
