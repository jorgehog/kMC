#pragma once


#include <sys/types.h>
#include <armadillo>

#include <libconfig_utils/libconfig_utils.h>

#include "site.h"


using namespace arma;


namespace kMC
{


class Reaction;

class KMCSolver
{
public:


    KMCSolver(const Setting & root);

    ~KMCSolver();


    void run();

    void initializeCrystal();

    void getRateVariables();

    uint getReactionChoice(double R);


    uint nNeighbors(uint & x, uint & y, uint & z)
    {
        return sites[x][y][z]->nNeighbors(0);
    }

    uint nNextNeighbors(uint & x, uint & y, uint & z)
    {
        return sites[x][y][z]->nNeighbors(1);
    }

    Site* getSite(const uint & i, const uint & j, const uint & k) const
    {
        return sites[i][j][k];
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

    const vector<double> & accuAllRates() const
    {
        return m_accuAllRates;
    }

    const vector<Reaction*> & allReactions() const
    {
        return m_allReactions;
    }

    const double & kTot() const
    {
        return m_kTot;
    }


private:

    double saturation;

    double RelativeSeedSize;


    Site**** sites;

    uint NX;
    uint NY;
    uint NZ;


    double m_kTot;
    vector<double> m_accuAllRates;

    vector<Reaction*> m_allReactions;

    double totalTime;


    uint m_nCycles;
    uint cycle;

    uint cyclesPerOutput;
    uint outputCounter;



    void initializeDiffusionReactions();

    void initializeSites();


    void dumpXYZ();

    void dumpOutput();

    static uint ptrCount;

};

}
