#include "sphericalEdge.h"

#include "../../kmcsolver.h"
#include "../../site.h"

using namespace kMC;

SphericalEdge::SphericalEdge(const uint dimension, const uint orientation) :
    Boundary(dimension, orientation, Boundary::SphericalEdge)
{

}

void SphericalEdge::update()
{
    for (Site* site : m_sitesBehindBoundary)
    {
        if (site->isActive())
        {
            solver()->despawnParticle(site->associatedParticle());
        }
    }
}

void SphericalEdge::initialize()
{
    m_sitesBehindBoundary.clear();
    m_interface.clear();

    uint x, y, z;

    double theta0, theta1, phi0, phi1;

    if (dimension() == 0)
    {
        theta0 = 0;
        theta1 = 1;

        if (orientation() == 0)
        {
            phi0 = 1;
            phi1 = 2;
        }
        else
        {
            phi0 = 0;
            phi1 = 1;
        }
    }

    else if (dimension() == 1)
    {
        theta0 = 0;
        theta1 = 1;

        if (orientation() == 0)
        {
            phi0 = 0.5;
            phi1 = 1.5;
        }
        else
        {
            phi0 = 1.5;
            phi1 = 2.5;
        }
    }

    else
    {
        phi0 = 0;
        phi1 = 2;

        if (orientation() == 0)
        {
            theta0 = 0.5;
            theta1 = 1;
        }
        else
        {
            theta0 = 0;
            theta1 = 0.5;
        }
    }

    theta0 *= datum::pi;
    theta1 *= datum::pi;
    phi0 *= datum::pi;
    phi1 *= datum::pi;

    vec thetas = linspace(theta0, theta1, 500);
    vec phis   = linspace(phi0, phi1, 500);

    vector<Site*> interfaceSites;
    for (const double &theta : thetas)
    {
        for (const double &phi : phis)
        {
            x = (NX() - 1.0)/2.0 + round(NX()*sin(theta)*cos(phi)/2.0);
            y = (NY() - 1.0)/2.0 + round(NY()*sin(theta)*sin(phi)/2.0);
            z = (NZ() - 1.0)/2.0 + round(NZ()*cos(theta)         /2.0);

            if (std::find(interfaceSites.begin(), interfaceSites.end(), solver()->getSite(x, y, z)) == interfaceSites.end())
            {
                m_interface.push_back({x, y, z});
                interfaceSites.push_back(solver()->getSite(x, y, z));
            }

        }
    }

    int sign = 1 - orientation()*2;

    for (uvec3 loc : m_interface)
    {
        while (loc(dimension()) != bound())
        {
            loc(dimension()) += sign;

            Site *site = solver()->getSite(loc(0), loc(1), loc(2));

            if (std::find(m_sitesBehindBoundary.begin(), m_sitesBehindBoundary.end(), site) == m_sitesBehindBoundary.end())
            {
                m_sitesBehindBoundary.push_back(site);
                if (!site->isActive())
                {
                    solver()->forceSpawnParticle(loc(0), loc(1), loc(2));
                }
            }

        }
    }

    solver()->dumpLAMMPS(1337);


}
