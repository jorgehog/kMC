#pragma once

#include "../periodic/periodic.h"


namespace kMC
{

class PeriodicShifted : public Periodic
{
public:
    PeriodicShifted(const uint dimension, const uint orientation, const int shift);

    const int &shift() const
    {
        return m_shift;
    }

    // Boundary interface
public:
    uint transformCoordinate(const int xi) const;
    int getDistanceBetween(int x1, int x2);

private:

    const int m_shift;
};

}
