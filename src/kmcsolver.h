
#pragma once


#include "site.h"

#include "debugger/debugger.h"

#include "RNG/kMCRNG.h"

#include <sys/types.h>
#include <armadillo>

#include <limits>

#include <libconfig_utils/libconfig_utils.h>


using namespace arma;


namespace kMC
{

class Reaction;

class KMCSolver
{
public:


    KMCSolver(const Setting & root);

    KMCSolver();

    ~KMCSolver();

    const static uint UNSET_UINT = std::numeric_limits<uint>::max();


    void mainloop();

    void reset();

    void initializeCrystal(const double relativeSeedSize);

    void initializeSolutionBath();

    void initializeSiteNeighborhoods()
    {
        forEachSiteDo([] (Site * site)
        {
            site->introduceNeighborhood();
        });
    }

    void initializeDiffusionReactions()
    {
        forEachSiteDo([] (Site * site)
        {
            site->initializeDiffusionReactions();
        });

    }

    void forEachSiteDo(function<void(Site * site)> applyFunction) const;

    void forEachSiteDo_sendIndices(function<void(Site *, uint, uint, uint)> applyFunction) const;


    void forEachActiveSiteDo(function<void(Site * site)> applyFunction) const;

    void forEachActiveSiteDo_sendIndices(function<void(Site *, uint, uint, uint)> applyFunction) const;


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

    const vector<Reaction*> & allPossibleReactions() const
    {
        return m_allPossibleReactions;
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


    void setTargetSaturation(const double saturation)
    {
        m_targetSaturation = saturation;
    }

    void setRNGSeed(uint seedState = Seed::fromTime, int defaultSeed = 0);



    void dumpXYZ();

    static void exit()
    {
        throw std::runtime_error("Exit failure.");
    }


    void clearSiteNeighborhoods()
    {
        forEachSiteDo([] (Site *site)
        {
            site->clearNeighborhood();
        });
    }

    void clearAllReactions()
    {
        forEachSiteDo([] (Site *site)
        {
            site->clearAllReactions();
        });
    }


    //OPTIMIZATION TMP

    void registerRateChange(const double &prevRate, const double &newRate)
    {


        if (prevRate == newRate)
        {
            return;
        }

        else if (prevRate == Reaction::UNSET_RATE)
        {

            KMCDebugger_Assert(newRate, !=, Reaction::UNSET_RATE);

            m_kTot += newRate;

        }

        else if (newRate == Reaction::UNSET_RATE)
        {
            m_kTot -= prevRate;
        }

        else
        {
            m_kTot += (newRate - prevRate);
        }

    }

    vector<Reaction*> prevUpdatedReacs;
    set<Reaction*> prevUpdatedReacsSet;


private:

    double m_targetSaturation = 0.01;

    Site**** sites;

    uint m_NX;
    uint m_NY;
    uint m_NZ;

    uvec3 m_N;


    double m_kTot;
    vector<double> m_accuAllRates;

    vector<Reaction*> m_allPossibleReactions;

    double totalTime;


    uint m_nCycles = 1000000;
    uint cycle;

    uint m_cyclesPerOutput;
    uint outputCounter;



    void initializeSites();

    void clearSites();

    void setBoxSize_KeepSites(const uvec3 &boxSizes);


    void dumpOutput();


    void checkRefCounter();

    void onConstruct();

    void finalizeObject();

    static uint refCounter;

    //OPTIMIZATION TMP
    vector<double> m_accuAllRates2;
    vector<Reaction*> m_allPossibleReactions2;


};

}
