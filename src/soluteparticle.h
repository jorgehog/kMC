#pragma once

#include "site.h"

#include <armadillo>

#include <vector>
#include <string>
#include <sstream>


namespace kMC
{

class Reaction;
class DiffusionReaction;

class SoluteParticle
{
public:

    SoluteParticle();

    ~SoluteParticle();

    void setSite(Site * site);

    void trySite(Site * site);

    void disableSite();

    void changeSite(Site *newSite);

    const Site *site() const
    {
        return m_site;
    }

    static void loadConfig(const Setting & setting);


    /*
     * Static non-trivial functions
     */


    static void selectUpdateFlags();

    static void popAffectedParticle(SoluteParticle *particle);

    static void updateAffectedParticles();




    /*
     *  Misc static property functions
     */


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
        return m_affectedParticles.find(this) != m_affectedParticles.end();
    }

    void markAsAffected();

    const vector<vector<SoluteParticle*> > & neighbouringParticles() const
    {
        return m_neighboringParticles;
    }


    const vector<SoluteParticle*> & neighbouringParticles(const uint level) const
    {
        return m_neighboringParticles.at(level);
    }

    const vector<SoluteParticle*> & closestNeighbors() const
    {
        return m_neighboringParticles.at(0);
    }



    const static set<SoluteParticle*, function<bool(SoluteParticle*, SoluteParticle*)> > & affectedParticles()
    {
        return m_affectedParticles;
    }

    const uint & x() const
    {
        return m_site->x();
    }

    const uint & y() const
    {
        return m_site->y();
    }

    const uint & z() const
    {
        return m_site->z();
    }

    const uint & r(const uint i) const
    {
        return m_site->r(i);
    }


    const vector<Reaction*> & reactions() const
    {
        return m_reactions;
    }


    DiffusionReaction* diffusionReactions(const uint i, const uint j, const uint k)
    {
        return m_diffusionReactions[i][j][k];
    }



    bool operator == (const SoluteParticle & other) const
    {
        return this == &other;
    }

    const string str() const
    {
        stringstream s;
        s << "SoluteParticle@(" << x() << ", " << y() << ", " << z() << ")";
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

    bool isNeighbor(SoluteParticle *particle, uint level);

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


    double potentialBetween(const SoluteParticle *other);


    inline void forEachNeighborDo(function<void (SoluteParticle*, const uint)> applyFunction) const;

    void forEachActiveReactionDo(function<void (Reaction*)> applyFunction) const;

    void forEachActiveReactionDo_sendIndex(function<void (Reaction*, uint)> applyFunction) const;


    const string info(int xr = 0, int yr = 0, int zr = 0, string desc = "X") const;

    void setZeroEnergy();

    void setVectorSizes();


private:

    static uvec4 m_totalParticles;

    static double m_totalEnergy;


    static set<SoluteParticle*, function<bool(SoluteParticle*, SoluteParticle*)> > m_affectedParticles;


    int m_particleState;


    vector<Reaction*> m_reactions;

    DiffusionReaction* m_diffusionReactions[3][3][3];


    Site *m_site;


    uvec m_nNeighbors;

    uint m_nNeighborsSum;

    double m_energy;


    vector<vector<SoluteParticle*> > m_neighboringParticles;


    void initializeDiffusionReactions();


    int detectParticleState();

    void setNewParticleState(int newState);

    void changeParticleState(int newState);


    const uint m_ID;

    static uint refCounter;


};


void SoluteParticle::forEachNeighborDo(function<void (SoluteParticle *, const uint)> applyFunction) const
{

    for (uint level = 0; level < Site::nNeighborsLimit(); ++level)
    {
        for (SoluteParticle *neighbor : m_neighboringParticles.at(level))
        {
            applyFunction(neighbor, level);
        }
    }

}



}
ostream& operator<<(ostream& os, const kMC::SoluteParticle& ss);

