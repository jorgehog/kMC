#include "interfacialstrain.h"

#include "../../debugger/debugger.h"

#include <algorithm>


using namespace kMC;

InterfacialStrain::InterfacialStrain(const Edge *interface,
                                     const double Es,
                                     const double r0,
                                     const uint nEdgeLayers) :
    Potential(),
    #ifndef KMC_NO_DEBUG
    m_trackedParticles(particleSet(SoluteParticle::compareFunc)),
    #endif
    m_Es(Es),
    m_r0(r0),
    m_interface(interface),
    m_nEdgeLayers(nEdgeLayers)
{
    if (m_nEdgeLayers == 0 || m_nEdgeLayers >= m_interface->span())
    {
        throw std::runtime_error("invalid number of interface layers.");
    }
}

InterfacialStrain::~InterfacialStrain()
{
    m_trackedParticles.clear();
    m_potential.clear();
}


void InterfacialStrain::initialize()
{
    m_potential.clear();

    m_potential.resize(2*m_interface->span() - 1);


    //first nEdgeLayer values are zero since there is no water between these sites.
    uint i = 0;
    for (double r = 0; r <= m_interface->span() - 1; r += 0.5, ++i)
    {
        m_potential.at(i) = strain(r - (int)m_nEdgeLayers);
    }

    KMCDebugger_Assert(i, ==, m_potential.size());

    if (m_interface->orientation() == 1)
    {
        reverse(m_potential.begin(), m_potential.end());
    }

}


double InterfacialStrain::valueAt(const double x, const double y, const double z)
{
    return strain(std::abs(selectXYZ(x, y, z) - (int)m_interface->bound()) - m_nEdgeLayers);
}

double InterfacialStrain::evaluateFor(SoluteParticle *particle)
{

    if (!isQualified(particle))
    {
        return 0;
    }

    return evaluateGivenQualified(particle);
}

double InterfacialStrain::evaluateSaddleFor(SoluteParticle *particle,
                                            const uint dx,
                                            const uint dy,
                                            const uint dz)
{
    const uint & r = particle->r(m_interface->dimension());
    int dr   = selectXYZ(dx, dy, dz) - 1;

    KMCDebugger_Assert((int)r + dr, >=, 0,"out of bounds.");
    KMCDebugger_Assert((int)r + dr, < , (int)m_interface->span(), "out of bounds.");

    return m_potential.at(2*r + dr);
}

double InterfacialStrain::onNeighborChange(SoluteParticle *particle,
                                           SoluteParticle *neighbor,
                                           const uint dx,
                                           const uint dy,
                                           const uint dz,
                                           int sign)
{

    (void) neighbor;
    (void) dx;
    (void) dy;
    (void) dz;
    (void) sign;


    if (!isQualified(particle))
    {

        if (isTracked(particle))
        {
            //not qualified and affected.
            m_trackedParticles.erase(particle);
        }

        return 0;
    }

    if (isTracked(particle))
    {
        return 0;
    }

    //qualified and not affected.
    m_trackedParticles.insert(particle);

    return evaluateGivenQualified(particle);

}

double InterfacialStrain::evaluateGivenQualified(SoluteParticle *particle)
{
    return m_potential.at(2*particle->r(m_interface->dimension()));
}

bool InterfacialStrain::isTracked(SoluteParticle *particle) const
{
    return m_trackedParticles.find(particle) != m_trackedParticles.end();
}

bool InterfacialStrain::isQualified(const SoluteParticle *particle) const
{
    return particle->isSurface();
}

double InterfacialStrain::strain(const double r) const
{
    //We demand at least one cell separation to induce a pressure
    if (r < 1)
    {
        return 0;
    }

    return m_Es*std::exp(-r/m_r0);
}

