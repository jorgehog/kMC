#include "inertwall.h"

#include "../../debugger/debugger.h"

#include "../../reactions/diffusion/diffusionreaction.h"

#include <algorithm>


using namespace kMC;

InertWall::InertWall(const Boundary *interface,
                     const double Es,
                     const double r0,
                     const double alpha_aw,
                     const double E_aw,
                     const double distanceFromEdge) :
    Potential(),
    m_Es(Es),
    m_r0(r0),
    m_alpha_aw(alpha_aw),
    m_E_aw(E_aw),
    m_interface(interface),
    m_distanceFromEdge(distanceFromEdge)
{
    if (!m_interface->type() == Boundary::Edge)
    {
        throw std::runtime_error("Invalid interface boundary type.");
    }

    if (!(E_aw/Es < 0))
    {
        throw std::runtime_error("E_aw and Es must have different signs");
    }
}

InertWall::~InertWall()
{
    m_trackedParticles.clear();
}


void InertWall::initialize()
{

}


double InertWall::valueAt(const double r, const double a, const double b)
{
    double rScaled = r + 1;

    if (m_interface->orientation() == 1)
    {
        rScaled = m_interface->span() - r - 1;
    }

    rScaled += m_distanceFromEdge;

    return strain(rScaled) + electroStatic(rScaled);
}

double InertWall::evaluateFor(SoluteParticle *particle)
{

    if (!isQualified(particle))
    {
        return 0;
    }

    //qualified. If already affected the set will not increase.
    m_trackedParticles.insert(particle->ID());

    return evaluateGivenQualified(particle);
}

double InertWall::evaluateSaddleFor(const DiffusionReaction *currentReaction)
{

    if (!isQualifiedSaddle(currentReaction))
    {
        return 0;
    }

    const uint & r = currentReaction->reactant()->r(m_interface->dimension());

    int dr = selectXYZ(currentReaction->path(0),
                       currentReaction->path(1),
                       currentReaction->path(2));

    KMCDebugger_Assert((int)r + dr, >=, 0,"out of bounds.");
    KMCDebugger_Assert((int)r + dr, <= , (int)m_interface->span(), "out of bounds.");

    return valueAt(r + dr/2, 0, 0);
}

double InertWall::onNeighborChange(SoluteParticle *particle,
                                   const SoluteParticle *neighbor,
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
            m_trackedParticles.erase(particle->ID());
            KMCDebugger_AssertBool(!isTracked(particle));

            return -evaluateGivenQualified(particle);
        }

        return 0;
    }

    if (isTracked(particle))
    {
        return 0;
    }

    //qualified. If already affected the set will not increase.
    m_trackedParticles.insert(particle->ID());

    return evaluateGivenQualified(particle);

}

double InertWall::evaluateGivenQualified(SoluteParticle *particle)
{
    return valueAt(particle->r(m_interface->dimension()), 0, 0);
}

bool InertWall::isTracked(SoluteParticle *particle) const
{
    return m_trackedParticles.find(particle->ID()) != m_trackedParticles.end();
}

bool InertWall::isQualified(const SoluteParticle *particle) const
{
    if (!particle->isSurface())
    {
        return false;
    }

    return hasCorrectOrientation(particle->x(), particle->y(), particle->z());
}

bool InertWall::isQualifiedSaddle(const DiffusionReaction *currentReaction) const
{
    const uint & r = currentReaction->reactant()->r(m_interface->dimension());

    if (m_interface->distanceFrom(r, true) <= (int)m_distanceFromEdge)
    {
        return false;
    }


    //Passing reaction site as phantom site to the calculations. This results in it being treated as invisible.
    return isQualified(currentReaction->reactant()) || hasCorrectOrientation(currentReaction->xD(),
                                                                             currentReaction->yD(),
                                                                             currentReaction->zD(),
                                                                             currentReaction->reactant()->site());
}

double InertWall::strain(const double r) const
{
    KMCDebugger_Assert(r, >, 0);

    return m_Es*std::exp(-r/m_r0)/r;
}

double InertWall::electroStatic(const double r) const
{
    KMCDebugger_Assert(r, >, 0);
    return m_E_aw/std::pow(r, m_alpha_aw);
}

bool InertWall::hasCorrectOrientation(const uint x, const uint y, const uint z, const Site *phantomSite) const
{
    int orientation = Site::detectSurfaceOrientation(x, y, z, m_interface->dimension(), phantomSite);

    if (orientation == 0)
    {
        return false;
    }

    //rescale orientations from -1 1 to 0 1
    return ((uint)(orientation + 1)/2 == m_interface->orientation());

}

