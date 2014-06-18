#pragma once

#include "../boundary.h"

namespace kMC
{


class ConcentrationWall : public Boundary
{
public:

    ConcentrationWall(const uint dimension, const uint orientation);

    ~ConcentrationWall();

    // Boundary interface
public:
    void update();
    void initialize() {}

private:


};

}
