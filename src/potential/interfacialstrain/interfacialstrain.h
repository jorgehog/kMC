#pragma once

#include "../potential.h"

#include <vector>

using namespace std;

namespace kMC
{

class Edge;

class InterfacialStrain : public Potential
{
public:
    InterfacialStrain(const Edge *interface, const double Es, const double r0, const uint nEdgeLayers);

    void initialize();

    double valueAt(const double x, const double y, const double z);

    double evaluateFor(SoluteParticle *particle);

    double evaluateSaddleFor(SoluteParticle *particle,
                             const uint dx,
                             const uint dy,
                             const uint dz);

    double onNeighborChange(SoluteParticle *neighbor,
                            const uint dx,
                            const uint dy,
                            const uint dz);

private:

    const double m_Es;

    const double m_r0;


    const Edge *m_interface;

    const uint m_nEdgeLayers;

    vector<double> m_potential;


    double strain(const double r) const;


};

}
