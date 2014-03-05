#pragma once


#include "../boundary.h"

namespace kMC
{

class Periodic : public Boundary
{
public:
    Periodic(uint orientation);

    // Boundary interface
public:
    uint transformCoordinate(const int xi)
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
