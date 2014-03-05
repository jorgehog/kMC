#include "periodic.h"


using namespace kMC;

Periodic::Periodic(uint orientation)
{
    switch (orientation) {
    case X:
        span = NX();
        break;
    case Y:
        span = NY();
        break;
    case Z:
        span = NZ();
        break;
    default:
        break;
    }

    delta.set_size(span);

    for(uint i = 0; i < span; ++i)
    {
        delta(i) = i;
        if (i > span/2)
        {
            delta(i) = -(int)(span - i);
        }
    }

}

Periodic::~Periodic()
{
    delta.clear();
}
