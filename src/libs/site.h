
#ifndef SITE_H
#define SITE_H

#include <vector>
#include <set>
#include <sys/types.h>
#include <armadillo>
#include <assert.h>

#include <libconfig_utils/libconfig_utils.h>

using namespace arma;

class KMCSolver;
class Reaction;

struct particleState {
    enum {
        crystal,
        solution,
        surface,
        any
    };
    const static vector<string> names;
};

class Site
{
public:

    Site(uint _x, uint _y, uint _z);

    ~Site();

    static const uint &nNeighborsLimit()
    {
        return m_nNeighborsLimit;
    }

    static const uint &neighborhoodLength()
    {
        return m_neighborhoodLength;
    }

    static const ucube &levelMatrix()
    {
        return m_levelMatrix;
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

    int getParticleState()
    {
        return m_particleState;
    }

    void setParticleState(int state);

    string getName() {

        string name;
        switch (m_particleState) {
        case particleState::crystal:
            name = "C";
            break;
        case particleState::solution:
            name = "P";
            break;
        case particleState::surface:
            name = "S";
            break;
        default:
            name = "X";
            break;
        }

        return name;
    }

    bool isLegalToSpawn();

    bool isCrystal() {
        return m_particleState == particleState::crystal;
    }

    bool isSurface() {
        return m_particleState == particleState::surface;
    }

    void crystallize();

    static void loadNeighborLimit(const Setting & setting);

    static uint getLevel(uint i, uint j, uint k);

    void addReaction(Reaction* reaction);

    void updateReactions();

    void spawnAsCrystal();

    void calculateRates();

    static void setSolverPtr(KMCSolver* solver);

    void distanceTo(const Site * other, int &dx, int &dy, int &dz, bool absolutes = false) const;


    uint nNeighbors(uint level = 0) const
    {
        return m_nNeighbors(level);
    }

    bool hasNeighboring(int state);

    void activate();

    void deactivate();

    const bool & active() const
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

    Site**** neighborHood() const
    {
        return m_neighborHood;
    }

    double energy() const
    {
        return m_energy;
    }

    void reset() {
        m_nNeighbors.zeros();
        m_totalEnergy -= m_energy;
        m_energy = 0;
    }

    static void resetAll() {
        m_totalActiveSites = 0;
        m_totalEnergy = 0;
        m_levelMatrix.reset();
        m_originTransformVector.reset();
    }

    void introduceNeighborhood();

    void propagateToNeighbors(int reqOldState, int newState);

    void informNeighborhoodOnChange(int change);

    void queueAffectedSites();



    void dumpInfo(int xr = 0, int yr = 0, int zr = 0);



    friend class testBed;

private:

    static uint m_nNeighborsLimit;
    static uint m_neighborhoodLength;

    static ucube m_levelMatrix;
    static ivec m_originTransformVector;

    static uint m_totalActiveSites;

    static uint NX;
    static uint NY;
    static uint NZ;

    static double m_totalEnergy;

    static set<Site*> affectedSites;
    static void updateAffectedSites();

    static KMCSolver* mainSolver;

    Site**** m_neighborHood;
    vector<Site*> m_allNeighbors;
    uvec m_nNeighbors;

    bool m_active;

    uint m_x;
    uint m_y;
    uint m_z;

    double m_energy;

    int m_particleState = particleState::solution;

    vector<Reaction*> m_siteReactions;

    vector<Reaction*> m_activeReactions;

};

#endif // SITE_H
