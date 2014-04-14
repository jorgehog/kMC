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

    uint transformCoordinate(const int xi) const;

    int getDistanceBetween(int x1, int x2);

    void initialize();

    void update() {}


private:

    ivec delta;

    void setupDelta();

};

}

