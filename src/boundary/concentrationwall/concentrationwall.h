#pragma once

#include "../boundary.h"

namespace kMC
{


class ConcentrationWall : public Boundary
{
public:

    ConcentrationWall(const uint dimension, const uint orientation);

    ~ConcentrationWall();

    static void setMaxEventsPrCycle(uint val)
    {
        m_maxEventsPrCycle = val;
    }

    static void setCooldown(uint val)
    {
        m_coolDown = val;
    }

    // Boundary interface
public:
    void update();
    void initialize() {}

private:

    //TMP
    static uint m_maxEventsPrCycle;
    static uint m_coolDown;
    static uint counter;

};

}
