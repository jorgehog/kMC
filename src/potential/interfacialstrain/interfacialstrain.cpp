#include "interfacialstrain.h"

#include "../../boundary/edge/edge.h"
#include "../../soluteparticle.h"
#include "../../debugger/debugger.h"

#include <algorithm>


using namespace kMC;

InterfacialStrain::InterfacialStrain(const Edge *interface,
                                     const double Es,
                                     const double r0,
                                     const uint nEdgeLayers) :
    Potential(),
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


void kMC::InterfacialStrain::initialize()
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


double kMC::InterfacialStrain::valueAt(const double x, const double y, const double z)
{
    double strainValue;

    switch (m_interface->dimension())
    {

    case Boundary::X:
        strainValue = strain(std::abs(x - (int)m_interface->bound()) - m_nEdgeLayers);
        break;

    case Boundary::Y:
        strainValue = strain(std::abs(y - (int)m_interface->bound()) - m_nEdgeLayers);
        break;

    case Boundary::Z:
        strainValue = strain(std::abs(z - (int)m_interface->bound()) - m_nEdgeLayers);
        break;

    default:
        break;
    }

    return strainValue;
}

double kMC::InterfacialStrain::evaluateFor(SoluteParticle *particle)
{
    return m_potential.at(2*particle->r(m_interface->dimension()));
}

double kMC::InterfacialStrain::evaluateSaddleFor(SoluteParticle *particle,
                                                 const uint dx,
                                                 const uint dy,
                                                 const uint dz)
{

}

double kMC::InterfacialStrain::onNeighborChange(SoluteParticle *neighbor,
                                                const uint dx,
                                                const uint dy,
                                                const uint dz)
{

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

