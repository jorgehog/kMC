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

    bool allowsTransitionTo(int state);

    bool isBlocked(){
        return nNeighbors(0) > 1;
    }

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

    void calculateRates();

    static void setSolverPtr(KMCSolver* solver);

    void distanceTo(const Site * other, int &dx, int &dy, int &dz, bool absolutes = false) const;


    uint nNeighbors(uint level = 0) const
    {
        return m_nNeighbors(level);
    }

    bool hasCrystalNeighbor();

    void activate()
    {

        assert(m_active == false && "activating active site.");
        assert(!isCrystal() && "A crystal should always be active");

        m_active = true;

        if (isSurface()) {
            cout << "I is surface" << endl;
            setParticleState(particleState::crystal);
        }

        updateReactions();
        calculateRates();

        informNeighborhoodOnChange(+1);

        m_totalActiveSites++;

    }

    void deactivate()
    {

        assert(m_active == true && "deactivating deactive site.");
        assert(!isSurface() && "A surface should never be active");

        m_active = false;

        //if we deactivate a crystal site, we have to potentially
        //reduce the surface by adding more sites as solution sites.
        //Site will change only if it is not surrounded by any crystals.
        if (isCrystal())
        {
            setParticleState(particleState::surface);
        }


        updateReactions();
        calculateRates();

        informNeighborhoodOnChange(-1);

        m_totalActiveSites--;

    }

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

    const std::vector<Reaction*> & activeReactions()
    {
        return m_activeReactions;
    }

    Site**** getNeighborhood() const
    {
        return neighborHood;
    }

    double getEnergy() const
    {
        return E;
    }

    void reset() {
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

    void countNeighbors();

    friend class testBed;


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

    int m_particleState   = particleState::solution;


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

    static KMCSolver* mainSolver;

    uvec m_nNeighbors;

    Site**** neighborHood;

    double E;

    uint m_x;
    uint m_y;
    uint m_z;

    bool m_active = false;

    std::vector<Reaction*> m_activeReactions;
    std::vector<Reaction*> m_siteReactions;

};

#endif // SITE_H
