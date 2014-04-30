#pragma once

#include "site.h"

#include <armadillo>

#include <vector>
#include <string>
#include <sstream>


#ifndef NDEBUG
#include <set>
#else
#include <unordered_set>
#endif

namespace kMC
{


#ifndef NDEBUG
typedef set<SoluteParticle*, function<bool(SoluteParticle*, SoluteParticle*)> > particleSet;
#else
typedef unordered_set<SoluteParticle*> particleSet;
#endif

class Reaction;
class DiffusionReaction;

class SoluteParticle
{
public:

    SoluteParticle();

    ~SoluteParticle();

    void setSite(const uint x, const uint y, const uint z);

    void trySite(const uint x, const uint y, const uint z);

    void resetSite();

    void disableSite();

    void changePosition(const uint x, const uint y, const uint z);

    const Site *site() const
    {
        return m_site;
    }

    static void loadConfig(const Setting & setting);


    /*
     * Static non-trivial functions
     */

    static void setMainSolver(KMCSolver* solver);


    static void selectUpdateFlags();

    static void popAffectedParticle(SoluteParticle *particle);

    static void updateAffectedParticles();




    /*
     *  Misc static property functions
     */


    static double getCurrentSolvantVolume();

    static double getCurrentConcentration();

    static double getCurrentRelativeCrystalOccupancy();


    /*
     * Init / Reset / clear static implementations
     */

    static void clearAll();

    static void clearAffectedParticles();

    static void setZeroTotalEnergy();


    /*
     * Misc. trivial functions
     */

    static KMCSolver * solver()
    {
        return m_solver;
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
        switch (i) {
        case 0:
            return m_x;
            break;
        case 1:
            return m_y;
            break;
        case 2:
            return m_z;
            break;
        default:
            break;
        }

        return m_x;

    }


    static const uint & nSurfaces()
    {
        return m_totalParticles.memptr()[ParticleStates::surface];
    }

    static uint nCrystals()
    {
        return m_totalParticles(ParticleStates::crystal);
    }

    static const uint & nSolutionParticles()
    {
        return m_totalParticles.memptr()[ParticleStates::solvant];
    }

    static uint nParticles()
    {
        return refCounter;
    }

    static const double & totalEnergy()
    {
        return m_totalEnergy;
    }

    static const uvec4 & totalParticlesVector()
    {
        return m_totalParticles;
    }

    static const uint & totalParticles(const uint i)
    {
        return m_totalParticles(i);
    }

    static bool isAffected(SoluteParticle *particle)
    {
        return m_affectedParticles.find(particle) != m_affectedParticles.end();
    }



    uint ID() const
    {
        return m_ID;
    }

    const int & particleState() const
    {
        return m_particleState;
    }


    uint nNeighbors(uint level = 0) const
    {
        return m_nNeighbors(level);
    }


    uint nNeighborsSum() const;


    double energy() const
    {
        return m_energy;
    }


    string particleStateName() const
    {
        return ParticleStates::names.at(m_particleState);
    }

    string particleStateShortName() const
    {
        return ParticleStates::shortNames.at(m_particleState);
    }

    uint nActiveReactions() const;

    bool isCrystal() const
    {
        return m_particleState == ParticleStates::crystal;
    }

    bool isSurface() const
    {
        return m_particleState == ParticleStates::surface;
    }

    bool isAffected()
    {
        return isAffected(this);
    }

    void markAsAffected();


    static const particleSet & affectedParticles()
    {
        return m_affectedParticles;
    }


    const vector<Reaction*> & reactions() const
    {
        return m_reactions;
    }


    DiffusionReaction* diffusionReactions(const uint i, const uint j, const uint k) const
    {
        return m_diffusionReactions[i][j][k];
    }



    bool operator == (const SoluteParticle & other) const
    {
        return this == &other;
    }

    bool operator < (const SoluteParticle * other) const
    {
        return this->ID() < other->ID();
    }

    bool operator > (const SoluteParticle * other) const
    {
        return this->ID() > other->ID();
    }

    const string str() const
    {
        stringstream s;
        s << "SoluteParticle" << m_ID << "@(" << x() << ", " << y() << ", " << z() << ")";
        return s.str();
    }

    const static uint & NX();

    const static uint & NY();

    const static uint & NZ();


    void clearAllReactions();



    /*
     * Non-trivial functions
     */

    void setParticleState(int newState);

    bool isLegalToSpawn() const;

    bool qualifiesAsCrystal() const
    {
        return nNeighbors() == Site::closestShellSize;
    }

    bool qualifiesAsSurface() const
    {
        return (nNeighbors() > 0) && !qualifiesAsCrystal();
    }

    bool isSolvant() const
    {
        return nNeighbors() == 0;
    }

    void reset();


    void addReaction(Reaction* reaction)
    {
        m_reactions.push_back(reaction);
    }

    void updateReactions();


    void setupAllNeighbors();

    void removeNeighbor(SoluteParticle *neighbor, uint level);

    void addNeighbor(SoluteParticle *neighbor, uint level);

    void _updateNeighborProps(const int sign, const SoluteParticle *neighbor, const uint level);


    void distanceTo(const SoluteParticle *other, int &dx, int &dy, int &dz, bool absolutes = false) const;

    double potentialBetween(const SoluteParticle *other);

    uint maxDistanceTo(const SoluteParticle *other) const;

    bool hasNeighboring(const int state) const
    {
        return Site::countNeighboring(m_x, m_y, m_z, state);
    }

    uint countNeighboring(const int state) const
    {
        return Site::countNeighboring(m_x, m_y, m_z, state);
    }


    void forEachActiveReactionDo(function<void (Reaction*)> applyFunction) const;

    void forEachActiveReactionDo_sendIndex(function<void (Reaction*, uint)> applyFunction) const;


    Site* neighborhood(const int dx, const int dy, const int dz) const
    {
        return Site::neighborhood(m_x, m_y, m_z, dx, dy, dz);
    }

    Site* neighborhood_fromIndex(const uint dx, const uint dy, const uint dz) const
    {
        return Site::neighborhood_fromIndex(m_x, m_y, m_z, dx, dy, dz);
    }

    void forEachNeighborSiteDo(function<void (Site *)> applyFunction)
    {
        Site::forEachNeighborDo(m_x, m_y, m_z, applyFunction);
    }

    void forEachNeighborSiteDo(function<void (Site *)> applyFunction) const
    {
        Site::forEachNeighborDo(m_x, m_y, m_z, applyFunction);
    }

    void forEachNeighborSiteDo_sendPath(function<void (Site *, int, int, int)> applyFunction)
    {
        Site::forEachNeighborDo_sendPath(m_x, m_y, m_z, applyFunction);
    }

    void forEachNeighborSiteDo_sendIndices(function<void (Site *, uint, uint, uint)> applyFunction)
    {
        Site::forEachNeighborDo_sendIndices(m_x, m_y, m_z, applyFunction);
    }


    const string info(int xr = 0, int yr = 0, int zr = 0, string desc = "X") const;

    void setZeroEnergy();

    void setVectorSizes();


private:

    static uvec4 m_totalParticles;

    static double m_totalEnergy;


    static particleSet m_affectedParticles;


    int m_particleState;


    vector<Reaction*> m_reactions;

    DiffusionReaction* m_diffusionReactions[3][3][3];


    Site *m_site;


    uint m_x;

    uint m_y;

    uint m_z;


    uvec m_nNeighbors;

    uint m_nNeighborsSum;

    double m_energy;


    void initializeDiffusionReactions();


    int detectParticleState();

    void setNewParticleState(int newState);

    void changeParticleState(int newState);


    const uint m_ID;

    static uint ID_count;

    static uint refCounter;

    static KMCSolver* m_solver;


};


}
ostream& operator<<(ostream& os, const kMC::SoluteParticle& ss);

