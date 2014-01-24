#ifndef SITE_H
#define SITE_H

#include <vector>
#include <sys/types.h>

class KMCSolver;
class Reaction;

class Site
{
public:
    Site(uint _x, uint _y, uint _z, KMCSolver* solver);

    void addReaction(Reaction* reaction);

    void updateReactions();

    void calculateRates();


    uint nNeighbours();

    uint nNextNeighbours();

    void activate() {
        m_active = true;
    }

    void deactivate() {
        m_active = false;
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

    std::vector<Reaction*> activeReactions() {
        return m_activeReactions;
    }

    //TEMPORARY SOLUTIONS
    double E;
    double En = 3;
    double Enn = 1;

private:

    KMCSolver* mainSolver;

    uint m_x;
    uint m_y;
    uint m_z;

    bool m_active = false;

    std::vector<Reaction*> m_activeReactions;
    std::vector<Reaction*> m_siteReactions;

};

#endif // SITE_H
