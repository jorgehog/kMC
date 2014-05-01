#pragma once


#include "../boundary.h"

namespace kMC
{

class Periodic : public Boundary
{
public:

    Periodic(const uint dimension, const uint orientation, const BoundaryTypes type = Boundary::Periodic);

    ~Periodic();


    // Boundary interface
public:

    virtual uint transformCoordinate(const int xi) const;

    virtual int getDistanceBetween(int x1, int x2);

    void initialize();

    void update() {}


private:

    ivec delta;

    void setupDelta();

};

}

