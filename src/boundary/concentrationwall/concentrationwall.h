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

    // Boundary interface
public:
    void update();
    void initialize() {}

private:

    //TMP
    static uint m_maxEventsPrCycle;
    static constexpr uint m_coolDown = 100u;
    static uint counter;

};

}
