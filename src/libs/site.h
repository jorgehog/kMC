#ifndef SITE_H
#define SITE_H

#include <vector>
#include <sys/types.h>
#include <armadillo>
#include <assert.h>

using namespace arma;

class KMCSolver;
class Reaction;

class Site
{
public:

    static const uint nNeighborsLimit = 2;
    static const uint neighborhoodLength;
    static ucube levelMatrix;

    static uint totalActiveSites;

    Site(uint _x, uint _y, uint _z, KMCSolver* solver);

    static uint getLevel(uint i, uint j, uint k);

    void addReaction(Reaction* reaction);

    void updateReactions();

    void calculateRates();


    uint nNeighbors(uint level = 0) {
        return m_nNeighbors(level);
    }

    void activate() {

        assert(m_active == false && "activating active site.");

        m_active = true;

        updateReactions();
        calculateRates();

        informNeighborhoodOnChange(+1);

        totalActiveSites++;

    }

    void deactivate() {

        assert(m_active == true && "deactivating deactive site.");

        m_active = false;

        updateReactions();
        calculateRates();

        informNeighborhoodOnChange(-1);

        totalActiveSites--;

    }

    bool active() {
        return m_active;
    }

    uint x()
    {
        return m_x;
    }

    uint y()
    {
        return m_y;
    }

    uint z()
    {
        return m_z;
    }

    const std::vector<Reaction*> & activeReactions() {
        return m_activeReactions;
    }

    Site**** getNeighborhood() {
        return neighborHood;
    }

    void introduceNeighborhood();

    void informNeighborhoodOnChange(int change);

    void countNeighbors();

    //TEMPORARY SOLUTIONS
    double E;
    double En = 3;
    double Enn = 1;

private:

    KMCSolver* mainSolver;

    uvec m_nNeighbors;

    Site**** neighborHood;

    uint m_x;
    uint m_y;
    uint m_z;

    bool m_active = false;

    std::vector<Reaction*> m_activeReactions;
    std::vector<Reaction*> m_siteReactions;


};

#endif // SITE_H
