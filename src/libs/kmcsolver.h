#pragma once


#include "site.h"

#include "RNG/kMCRNG.h"

#include <sys/types.h>
#include <armadillo>

#include <climits>

#include <libconfig_utils/libconfig_utils.h>


using namespace arma;


namespace kMC
{

class Reaction;

class KMCSolver
{
public:


    KMCSolver(const Setting & root);

    ~KMCSolver();

    const static uint UNSET_UINT = (uint)ULLONG_MAX;


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

    Site* getSite(const uint i, const uint j, const uint k) const
    {
        return sites[i][j][k];
    }

    const uint &NX () const
    {
        return m_NX;
    }

    const uint &NY () const
    {
        return m_NY;
    }

    const uint &NZ () const
    {
        return m_NZ;
    }

    const uint &N(const uint i) const
    {
        return m_N(i);
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


    //Set functions

    void setBoxSize(const uvec3 boxSize);

    void setNumberOfCycles(const uint nCycles)
    {
        m_nCycles = nCycles;
    }

    void setCyclesPerOutput(const uint cyclesPerOutput)
    {
        m_cyclesPerOutput = cyclesPerOutput;
    }

    void setRelativeSeedSize(const double relativeSeedSize)
    {

        if (relativeSeedSize > 1.0)
        {
            cerr << "The seed size cannot exceed the box size." << endl;
            exit(1);
        }

        m_relativeSeedSize = relativeSeedSize;


    }

    void setSaturation(const double saturation)
    {
        m_saturation = saturation;
    }

    void setRNGSeed(uint seedState, int defaultSeed);



    void dumpXYZ();


private:

    double m_saturation;

    double m_relativeSeedSize;


    Site**** sites;

    uint m_NX;
    uint m_NY;
    uint m_NZ;

    uvec3 m_N;


    double m_kTot;
    vector<double> m_accuAllRates;

    vector<Reaction*> m_allReactions;

    double totalTime;


    uint m_nCycles;
    uint cycle;

    uint m_cyclesPerOutput;
    uint outputCounter;



    void initializeDiffusionReactions();

    void initializeSites();

    void initializeSiteNeighborhood();

    void clearSites();


    void dumpOutput();


    static uint ptrCount;

};

}
