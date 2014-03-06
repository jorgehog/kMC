#pragma once


#include "../boundary.h"

namespace kMC
{

class Periodic : public Boundary
{
public:

    Periodic(const uint dimension, const uint orientation);

    ~Periodic();


    // Boundary interface
public:

    uint transformCoordinate(const int xi) const
    {
        return (xi + span())%span();
    }

    int getDistanceBetween(int x1, int x2)
    {
        return delta(transformCoordinate(Boundary::getDistanceBetween(x1, x2)));
    }

    void initialize() {}


private:

    ivec delta;


};

}

