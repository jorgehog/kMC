#pragma once

#include "../potential.h"

#include "../../boundary/edge/edge.h"

#include "../../soluteparticle.h"

using namespace std;

namespace kMC
{

class StressedSurface : public Potential
{
public:

    StressedSurface(const Edge *interface, const double Es, const double r0, const uint nEdgeLayers);

    ~StressedSurface();


    void initialize();

    double valueAt(const double x, const double y, const double z);

    double evaluateFor(SoluteParticle *particle);

    double evaluateSaddleFor(SoluteParticle *particle,
                             const uint dx,
                             const uint dy,
                             const uint dz);

    double onNeighborChange(SoluteParticle *particle,
                            const SoluteParticle *neighbor,
                            const uint dx,
                            const uint dy,
                            const uint dz,
                            int sign);

    bool isTracked(SoluteParticle *particle) const;

    bool isQualified(const SoluteParticle *particle) const;

private:

    particleSet m_trackedParticles;

    const double m_Es;

    const double m_r0;


    const Edge *m_interface;

    const uint m_nEdgeLayers;

    vector<double> m_potential;

    double evaluateGivenQualified(SoluteParticle *particle);


    double strain(const double r) const;

    template<typename T>
    T selectXYZ(const T &x, const T &y, const T &z) const
    {
        switch (m_interface->dimension())
        {

        case 0:
            return x;
            break;

        case 1:
            return y;
            break;

        case 2:
            return z;
            break;

        default:
            break;
        }

        return x;
    }


};

}
