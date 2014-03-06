#pragma once

#include <vector>
#include <set>
#include <sys/types.h>
#include <armadillo>
#include <assert.h>

#include <libconfig_utils/libconfig_utils.h>

using namespace arma;


namespace kMC
{


class KMCSolver;
class Reaction;
class DiffusionReaction;
class Boundary;

struct ParticleStates
{
    enum
    {
        crystal,
        solution,
        surface
    };

    const static vector<string> names;
    const static vector<string> shortNames;

};

class Site
{
public:

    Site(uint _x, uint _y, uint _z);

    ~Site();

    /*
     * Static non-trivial functions
     */

    static void setMainSolver(KMCSolver* solver);

    static void loadConfig(const Setting & setting);

    static void initializeBoundaries();

    static void updateBoundaries();

    static void updateAffectedSites();

    static void selectUpdateFlags();

    static uint findLevel(uint i, uint j, uint k);

    static void resetAll();

    /*
     * Non-trivial functions
     */

    void setParticleState(int state);


    bool isLegalToSpawn();

    void spawnAsFixedCrystal();


    void activate();

    void deactivate();


    void addReaction(Reaction* reaction);

    void updateReactions();

    void calculateRates();


    void introduceNeighborhood();

    bool hasNeighboring(int state) const;

    void propagateToNeighbors(int reqOldState, int newState);

    void informNeighborhoodOnChange(int change);


    void distanceTo(const Site * other, int &dx, int &dy, int &dz, bool absolutes = false) const;

    uint maxDistanceTo(const Site * other) const;

    double potentialBetween(const Site * other);

    void setDirectUpdateFlags();

    void queueAffectedSites();

    void setZeroEnergy();

    void reset()
    {
        m_nNeighbors.zeros();
        m_totalEnergy -= m_energy;
        m_energy = 0;
        m_nNeighborsSum = 0;
    }


    const string info(int xr = 0, int yr = 0, int zr = 0, string desc = "X") const;


    /*
     * Misc. trivial functions
     */


    static const uint &nNeighborsLimit()
    {
        return m_nNeighborsLimit;
    }

    static const uint &neighborhoodLength()
    {
        return m_neighborhoodLength;
    }

    static const uint &levelMatrix(const uint i, const uint j, const uint k)
    {
        return m_levelMatrix(i, j, k);
    }

    static const ivec &originTransformVector()
    {
        return m_originTransformVector;
    }

    static const uint & totalActiveSites()
    {
        return m_totalActiveSites;
    }

    static const double & totalEnergy()
    {
        return m_totalEnergy;
    }

    static const Boundary * boundaries(const uint xyz, const uint loc)
    {
        return m_boundaries(xyz, loc);
    }


    const int & particleState() const
    {
        return m_particleState;
    }

    string particleStateName() const
    {
        return ParticleStates::names.at(m_particleState);
    }

    string particleStateShortName() const
    {
        return ParticleStates::shortNames.at(m_particleState);
    }

    uint nNeighbors(uint level = 0) const
    {
        return m_nNeighbors(level);
    }

    uint nNeighborsSum() const;

    bool isCrystal() const
    {
        return m_particleState == ParticleStates::crystal;
    }

    bool isSurface() const
    {
        return m_particleState == ParticleStates::surface;
    }

    const bool & isActive() const
    {
        return m_active;
    }


    const uint & x() const
    {
        return m_x;
    }

    const uint & y() const
    {
        return m_y;
    }

    const uint & z() const
    {
        return m_z;
    }


    const vector<Reaction*> & activeReactions() const
    {
        return m_activeReactions;
    }

    const vector<Reaction*> & siteReactions() const
    {
        return m_siteReactions;
    }

    const vector<Site*> & allNeighbors() const
    {
        return m_allNeighbors;
    }

    const static set<Site*> & affectedSites()
    {
        return m_affectedSites;
    }

    Site* neighborHood(const uint x, const uint y, const uint z) const
    {
        return m_neighborHood[x][y][z];
    }

    double energy() const
    {
        return m_energy;
    }

    const bool & isFixedCrystalSeed()
    {
        return m_isFixedCrystalSeed;
    }

    bool operator == (const Site & other) const
    {
        return this == &other;
    }

    const string str() const
    {
        stringstream s;
        s << "Site@(" << x() << ", " << y() << ", " << z() << ")";
        return s.str();
    }



private:

    static field<Boundary*> m_boundaries;

    static uint m_nNeighborsLimit;
    static uint m_neighborhoodLength;

    static ucube m_levelMatrix;
    static ivec m_originTransformVector;

    static uint m_totalActiveSites;

    static uint NX;
    static uint NY;
    static uint NZ;

    static double m_totalEnergy;

    static set<Site*> m_affectedSites;

    static KMCSolver* mainSolver;

    Site**** m_neighborHood;
    vector<Site*> m_allNeighbors;
    uvec m_nNeighbors;
    uint m_nNeighborsSum;

    bool m_active;

    bool m_isFixedCrystalSeed;

    uint m_x;
    uint m_y;
    uint m_z;

    double m_energy;

    int m_particleState = ParticleStates::solution;

    vector<Reaction*> m_siteReactions;

    vector<Reaction*> m_activeReactions;


};

}

ostream& operator<<(ostream& os, const kMC::Site& ss);
