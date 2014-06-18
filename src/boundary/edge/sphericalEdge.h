#include "../boundary.h"

namespace kMC
{


class SphericalEdge : public Boundary
{
public:

    SphericalEdge(const uint dimension, const uint orientation);

    void update();
    void initialize();

    uint interfaceValue(const uint w, const uint l) const
    {
        return m_interface(w, l);
    }

private:

    umat m_interface;

    void removeOutsiders();

};


}
