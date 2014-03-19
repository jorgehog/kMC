#pragma once

#include "../boundary.h"

namespace kMC
{


class ConcentrationWall : public Boundary
{
public:

    ConcentrationWall(const uint dimension, const uint orientation);

    ~ConcentrationWall();

    static void setMinDistanceFromSite(const uint minDistanceFromSite)
    {
        m_minDistanceFromSurface = minDistanceFromSite;
    }

    static void setMaxEventsPrCycle(uint val)
    {
        m_maxEventsPrCycle = val;
    }

    // Boundary interface
public:
    void update();
    void initialize();
    void finalize();

private:

    static uint m_minDistanceFromSurface;

    static uint m_maxEventsPrCycle;

};

}
