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

private:

    uint span;

    enum Orientations
    {
        X,
        Y,
        Z
    };

};

}
