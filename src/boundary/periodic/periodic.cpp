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
    return Site::nNeighborsLimit() + (xi - Site::nNeighborsLimit() + span())%span();
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
