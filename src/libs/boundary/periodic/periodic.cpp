#include "periodic.h"


using namespace kMC;

Periodic::Periodic(const uint dimension, const uint orientation) :
    Boundary(dimension, orientation)
{

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

Periodic::~Periodic()
{
    delta.clear();
}