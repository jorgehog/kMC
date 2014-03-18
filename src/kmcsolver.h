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

    void reset();

    void initializeCrystal();

    void initializeSiteNeighborhoods();

    void initializeDiffusionReactions();


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

    const uvec3 NVec() const
    {
        return m_N;
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

    const double & targetSaturation()
    {
        return m_targetSaturation;
    }


    //Set functions

    void setBoxSize(const uvec3 boxSize, bool check = true, bool keepSystem = false);

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
            KMCSolver::exit();
        }

        m_relativeSeedSize = relativeSeedSize;


    }

    void setTargetSaturation(const double saturation)
    {
        m_targetSaturation = saturation;
    }

    void setRNGSeed(uint seedState, int defaultSeed);



    void dumpXYZ();

    static void exit()
    {
        throw std::runtime_error("Exit failure.");
    }

    //SHOULD BE IN SITE

    void clearSiteNeighborhoods();

    void clearAllReactions();


private:

    double m_targetSaturation;

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



    void initializeSites();

    void clearSites();

    void setBoxSize_KeepSites(const uvec3 &boxSizes);


    void dumpOutput();

    void checkRefCounter();

    static uint refCounter;


};

}