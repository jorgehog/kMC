#pragma once


#include "../boundary.h"

namespace kMC
{

class Periodic : public Boundary
{
public:

    Periodic(uint orientation);

    ~Periodic();


    // Boundary interface
public:

    uint transformCoordinate(const int xi) const
    {
        return (xi + span)%span;
    }

    int getDistanceBetween(int x1, int x2)
    {
        return delta(transformCoordinate(Boundary::getDistanceBetween(x1, x2)));
    }

private:

    uint span;

    ivec delta;

    enum Orientations
    {
        X,
        Y,
        Z
    };


};

}

