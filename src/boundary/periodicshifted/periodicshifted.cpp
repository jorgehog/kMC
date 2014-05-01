#include "periodicshifted.h"

using namespace kMC;


PeriodicShifted::PeriodicShifted(const uint dimension, const uint orientation, const int shift) :
    Periodic(dimension, orientation, Boundary::PeriodicShifted),
    m_shift(shift)
{

}

uint PeriodicShifted::transformCoordinate(const int xi) const
{
    uint xb = Boundary::transformCoordinate(xi);

    if (isBlocked(xb))
    {
        return Periodic::transformCoordinate(xi) + m_shift;
    }

    return xb;
}

int PeriodicShifted::getDistanceBetween(int x1, int x2)
{
    return Periodic::getDistanceBetween(x1, x2);
}
