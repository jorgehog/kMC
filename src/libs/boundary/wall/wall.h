#pragma once

#include "../boundary.h"

namespace kMC
{

class Wall : public Boundary
{
public:
    Wall(const uint dimension, const uint orientation);

    ~Wall();

    // Boundary interface
public:
    void update() {}
    void initialize();
};

}
