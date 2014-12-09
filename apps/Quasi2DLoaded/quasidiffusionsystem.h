#pragma once

#include <kMC>

#include "movingwall.h"

using namespace arma;

namespace kMC
{

class QuasiDiffusionSystem
{
private:

    ivec &m_heights;
    const double m_Eb;

    MovingWall &m_wallEvent;

public:

    QuasiDiffusionSystem(ivec &heights, const double Eb, MovingWall &wallEvent) :
        m_heights(heights),
        m_Eb(Eb),
        m_wallEvent(wallEvent)
    {

    }

    void setHeights(const ivec &heights);

    const double &Eb() const
    {
        return m_Eb;
    }

    const MovingWall &wallEvent() const
    {
        return m_wallEvent;
    }

    uint leftSite(const uint site, const uint n) const
    {
        return (site + size() - n)%size();
    }

    uint rightSite(const uint site, const uint n) const
    {
        return (site + n)%size();
    }

    bool connectedLeft(const uint leftSite, const int myHeight) const
    {
        return m_heights(leftSite) >= myHeight;
    }

    bool connectedRight(const uint rightSite, const int myHeight) const
    {
        return m_heights(rightSite) >= myHeight;
    }

    uint nNeighbors(const uint leftsite, const uint rightsite, const uint centersite) const
    {
        bool leftHug = connectedLeft(leftsite, m_heights(centersite));
        bool rightHug = connectedRight(rightsite, m_heights(centersite));

        if (leftHug && rightHug)
        {
            return 3;
        }

        else if (leftHug || rightHug)
        {
            return 2;
        }

        else
        {
            return 1;
        }
    }

    uint nNeighbors(const uint site) const
    {
        return nNeighbors(leftSite(site, 1), rightSite(site, 1), site);
    }

    ivec heights() const
    {
        return m_heights;
    }

    const int &heights(const uint i) const
    {
        return m_heights(i);
    }


    void registerHeightChange(const uint site, const int change)
    {
        m_heights(site) += change;
    }

    const uint &size() const
    {
        return m_heights.n_elem;
    }

    void markAsAffected(SoluteParticle *particle)
    {
        m_wallEvent.markAsAffected(particle);
    }

};

}
