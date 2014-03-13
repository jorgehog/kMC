#pragma once

#include "../boundary.h"

namespace kMC
{

class Edge : public Boundary
{
public:
    Edge(const uint dimension, const uint orientation) :
        Boundary(dimension, orientation, Boundary::Edge)
    {

    }

    // Boundary interface
public:
    void update() {}
    void initialize() {}
};

}
