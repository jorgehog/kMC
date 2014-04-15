#include "periodic.h"

#include "../../site.h"

using namespace kMC;

Periodic::Periodic(const uint dimension, const uint orientation) :
    Boundary(dimension, orientation, Boundary::Periodic)
{

}

Periodic::~Periodic()
{
    delta.clear();
}

uint Periodic::transformCoordinate(const int xi) const
{
    return (xi + span())%span();
}

int Periodic::getDistanceBetween(int x1, int x2)
{
    return delta(transformCoordinate(Boundary::getDistanceBetween(x1, x2)));
}

void Periodic::initialize()
{
    setupDelta();
}

void Periodic::setupDelta()
{
    delta.clear();
    delta.set_size(span());

    for(uint i = 0; i < span(); ++i)
    {
        delta(i) = i;
        if (i > span()/2)
        {
            delta(i) = -(int)(span() - i);
        }
    }
}
