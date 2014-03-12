#pragma once

#include "../boundary.h"

namespace kMC
{

class Surface : public Boundary
{
public:
    Surface(const uint dimension, const uint orientation);

    ~Surface();

    // Boundary interface
public:
    void update() {}
    void initialize();
};

}
