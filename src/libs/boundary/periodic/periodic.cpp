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
}

Periodic::~Periodic()
{

}
