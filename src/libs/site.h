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


    static uint m_nNeighborsLimit;
    static uint m_neighborhoodLength;

    static ucube m_levelMatrix;
    static ivec m_originTransformVector;

    static uint totalActiveSites;

    Site(uint _x, uint _y, uint _z);

    static void loadNeighborLimit(const Setting & setting);

    static uint getLevel(uint i, uint j, uint k);

    void addReaction(Reaction* reaction);

    void updateReactions();

    void calculateRates();

    static void setSolverPtr(KMCSolver* solver);


    uint nNeighbors(uint level = 0)
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

        totalActiveSites++;

    }

    void deactivate()
    {

        assert(m_active == true && "deactivating deactive site.");

        m_active = false;

        updateReactions();
        calculateRates();

        informNeighborhoodOnChange(-1);

        totalActiveSites--;

    }

    void updateEnergy(Site *changedSite, int change);

    const bool & active()
    {
        return m_active;
    }

    const uint & x()
    {
        return m_x;
    }

    const uint & y()
    {
        return m_y;
    }

    const uint & z()
    {
        return m_z;
    }

    const std::vector<Reaction*> & activeReactions()
    {
        return m_activeReactions;
    }

    Site**** getNeighborhood()
    {
        return neighborHood;
    }

    double getEnergy()
    {
        return E;
    }

    void introduceNeighborhood();

    void informNeighborhoodOnChange(int change);

    void countNeighbors();

    friend class testBed;


    static const uint &nNeighborsLimit();

    static const uint &neighborhoodLength();

    static const ucube &levelMatrix();

    static const ivec &originTransformVector();


private:

    static uint NX;
    static uint NY;
    static uint NZ;

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
