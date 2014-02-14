
#ifndef SITE_H
#define SITE_H

#include <vector>
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

    int getParticleState() {
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

    void linkSiteDependency(const Site * site) {
        m_siteDependencies.push_back(site);
    }

    void crystallize();

    static void loadNeighborLimit(const Setting & setting);

    static uint getLevel(uint i, uint j, uint k);

    void addReaction(Reaction* reaction);

    void updateReactions();

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


    void updateEnergy(Site *changedSite, int change);


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

    const vector<const Site*> & siteDependencies()  const
    {
        return m_siteDependencies;
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

    double getEnergy() const
    {
        return E;
    }

    void reset() {
        m_siteDependencies.clear();
        m_nNeighbors.zeros();
        m_totalEnergy -= E;
        E = 0;
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

    void updateNeighborReactions();

    static const uint &nNeighborsLimit();

    static const uint &neighborhoodLength();

    static const ucube &levelMatrix();

    static const ivec &originTransformVector();

    static const uint & totalActiveSites() {
        return m_totalActiveSites;
    }

    static const double & totalEnergy() {
        return m_totalEnergy;
    }

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

    static vector<Site*> affectedSites;

    static KMCSolver* mainSolver;

    uvec m_nNeighbors;

    Site**** m_neighborHood;

    double E;

    uint m_x;
    uint m_y;
    uint m_z;

    bool m_active = false;

    int m_particleState = particleState::solution;

    vector<Reaction*> m_activeReactions;
    vector<Reaction*> m_siteReactions;

    vector<Site*> m_allNeighbors;
    vector<const Site*> m_siteDependencies;

};

#endif // SITE_H
