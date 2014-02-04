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

class Site
{
public:



    Site(uint _x, uint _y, uint _z);

    ~Site();

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

    void activate()
    {

        assert(m_active == false && "activating active site.");

        m_active = true;

        updateReactions();
        calculateRates();

        informNeighborhoodOnChange(+1);

        m_totalActiveSites++;

    }

    void deactivate()
    {

        assert(m_active == true && "deactivating deactive site.");

        m_active = false;

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
