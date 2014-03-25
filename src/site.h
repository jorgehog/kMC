#pragma once

#include "particlestates.h"

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


    /*
     *  Misc static property functions
     */

    static uint getLevel(uint i, uint j, uint k);

    static double getCurrentSolutionDensity();

    static double getCurrentRelativeCrystalOccupancy();

    static umat getCurrentCrystalBoxTopology();



    /*
     * Init / Reset / clear static implementations
     */

    static void setInitialNNeighborsLimit(const uint & nNeighborsLimit, bool check = true);

    static void setInitialBoundaries(const umat & boundaryMatrix);

    static void setInitialBoundaries(const int boundaryType);

    static void setInitialNNeighborsToCrystallize(const uint & nNeighborsToCrystallize);


    static void resetBoundariesTo(const umat & boundaryMatrix);

    static void resetBoundariesTo(const int boundaryType);

    static void resetNNeighborsLimitTo(const uint & nNeighborsLimit, bool check = true);

    static void resetNNeighborsToCrystallizeTo(const uint & nNeighborsToCrystallize);


    static void clearAll();

    static void clearAffectedSites();

    static void clearBoundaries();

    static void finalizeBoundaries();


    static void setZeroTotalEnergy();

    /*
     * Non-trivial functions
     */


    void setParticleState(int newState);

    bool isLegalToSpawn();


    bool qualifiesAsCrystal();

    bool qualifiesAsSurface();


    void spawnAsFixedCrystal();

    void spawnAsCrystal();


    void blockCrystallizationOnSite();

    void allowCrystallizationOnSite();


    void crystallize();

    void decrystallize();


    void activate();

    void deactivate();

    void flipActive();

    void flipDeactive();


    void addReaction(Reaction* reaction)
    {
        m_reactions.push_back(reaction);
    }

    void calculateRates();


    void initializeDiffusionReactions();

    void introduceNeighborhood();


    bool hasNeighboring(int state, int range) const;

    uint countNeighboring(int state, int range) const;


    void propagateToNeighbors(int reqOldState, int newState, int range);

    void informNeighborhoodOnChange(int change);


    void distanceTo(const Site * other, int &dx, int &dy, int &dz, bool absolutes = false) const;

    uint maxDistanceTo(const Site * other) const;

    double potentialBetween(const Site * other);

    void setNeighboringDirectUpdateFlags();

    void queueAffectedSites();

    void setZeroEnergy();

    void reset();

    void clearNeighborhood();


    const string info(int xr = 0, int yr = 0, int zr = 0, string desc = "X") const;


    void forEachNeighborDo(function<void (Site *)> applyFunction) const;

    void forEachNeighborDo_sendIndices(function<void (Site *, uint, uint, uint)> applyFunction) const;


    void forEachActiveReactionDo(function<void (Reaction*)> applyFunction) const;

    void forEachActiveReactionDo_sendIndex(function<void (Reaction*, uint)> applyFunction) const;


    /*
     * Misc. trivial functions
     */


    static const uint & nSurfaces()
    {
        return m_totalDeactiveParticles.memptr()[ParticleStates::surface];
    }

    static uint nCrystals()
    {
        return m_totalActiveParticles(ParticleStates::crystal) + m_totalActiveParticles(ParticleStates::fixedCrystal);
    }

    static const uint & nSolutionParticles()
    {
        return m_totalActiveParticles.memptr()[ParticleStates::solution];
    }


    static const uint &boundaryTypes(const uint i, const uint j = 0)
    {
        return m_boundaryTypes(i, j);
    }

    static const uint &nNeighborsToCrystallize()
    {
        return m_nNeighborsToCrystallize;
    }

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

    static uint originTransformVector(const uint i)
    {
        return m_originTransformVector(i);
    }

    static const uint & totalActiveSites()
    {
        return m_totalActiveSites;
    }

    static const uvec4 & totalActiveParticlesVector()
    {
        return m_totalActiveParticles;
    }

    static const uvec4 & totalDeactiveParticlesVector()
    {
        return m_totalDeactiveParticles;
    }

    static const uint & totalActiveParticles(const uint i)
    {
        return m_totalActiveParticles(i);
    }

    static const uint & totalDeactiveParticles(const uint i)
    {
        return m_totalDeactiveParticles(i);
    }

    static const double & totalEnergy()
    {
        return m_totalEnergy;
    }

    static const Boundary * boundaries(const uint xyz, const uint loc)
    {
        return m_boundaries(xyz, loc);
    }

    static const field<Boundary*> & boundaryField()
    {
        return m_boundaries;
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

    uint nActiveReactions() const;

    uint nNeighborsSum() const;

    bool isCrystal() const
    {
        return ((m_particleState == ParticleStates::crystal) || (m_particleState == ParticleStates::fixedCrystal));
    }

    bool isSurface() const
    {
        return m_particleState == ParticleStates::surface;
    }

    const bool & isActive() const
    {
        return m_active;
    }

    bool isAffected()
    {
        return m_affectedSites.find(this) != m_affectedSites.end();
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

    const uint & r(const uint i) const
    {
        return m_r(i);
    }

    const vector<Reaction*> & reactions() const
    {
        return m_reactions;
    }

    const static set<Site*> & affectedSites()
    {
        return m_affectedSites;
    }

    Site* neighborhood(const uint x, const uint y, const uint z) const
    {
        return m_neighborhood[x][y][z];
    }

    double energy() const
    {
        return m_energy;
    }

    const bool & isFixedCrystalSeed()
    {
        return m_isFixedCrystalSeed;
    }

    const bool & cannotCrystallize()
    {
        return m_cannotCrystallize;
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

    const static uint & NX();

    const static uint & NY();

    const static uint & NZ();

    //should be in site.. all sites in sites and jazz jazz..

    void clearAllReactions();


private:

    static field<Boundary*> m_boundaries;

    static field<const Setting*> m_boundaryConfigs;

    static umat m_boundaryTypes;


    static uint m_nNeighborsLimit;

    static uint m_neighborhoodLength;


    static uint m_nNeighborsToCrystallize;


    static ucube m_levelMatrix;

    static ivec m_originTransformVector;


    static uint m_totalActiveSites;

    static uvec4 m_totalActiveParticles;

    static uvec4 m_totalDeactiveParticles;

    static double m_totalEnergy;


    static set<Site*> m_affectedSites;

    static KMCSolver* m_solver;


    Site**** m_neighborhood;

    uvec m_nNeighbors;

    uint m_nNeighborsSum;


    bool m_active;

    bool m_isFixedCrystalSeed;

    bool m_cannotCrystallize;


    const uint m_x;

    const uint m_y;

    const uint m_z;

    const uvec3 m_r;


    double m_energy;


    int m_particleState = ParticleStates::solution;


    vector<Reaction*> m_reactions;

    DiffusionReaction* m_diffusionReactions[3][3][3];


    void setNewParticleState(int newState);

    void deactivateFixedCrystal();


};

}

ostream& operator<<(ostream& os, const kMC::Site& ss);
