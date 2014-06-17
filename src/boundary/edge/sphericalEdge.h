#include "../boundary.h"

namespace kMC
{


class SphericalEdge : public Boundary
{
public:

    SphericalEdge(const uint dimension, const uint orientation);

    // Boundary interface
public:
    void update();
    void initialize();

private:

    umat m_interface;

};


}
