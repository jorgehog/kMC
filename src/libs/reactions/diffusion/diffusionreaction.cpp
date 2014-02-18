#include "diffusionreaction.h"
#include "../../kmcsolver.h"

DiffusionReaction::DiffusionReaction(Site *destination) :
    Reaction(),
    destination(destination)
{

}


void DiffusionReaction::loadConfig(const Setting &setting)
{

    rPower = getSurfaceSetting<double>(setting, "rPower");
    scale  = getSurfaceSetting<double>(setting, "scale");

    m_potential.set_size(Site::neighborhoodLength(),
                         Site::neighborhoodLength(),
                         Site::neighborhoodLength());

    for (uint i = 0; i < Site::neighborhoodLength(); ++i)
    {
        for (uint j = 0; j < Site::neighborhoodLength(); ++j)
        {
            for (uint k = 0; k < Site::neighborhoodLength(); ++k)
            {

                if (i == Site::nNeighborsLimit() && j == Site::nNeighborsLimit() && k == Site::nNeighborsLimit())
                {
                    m_potential(i, j, k) = 0;
                    continue;
                }

                m_potential(i, j, k) = 1.0/pow(pow(Site::originTransformVector()(i), 2)
                                               + pow(Site::originTransformVector()(j), 2)
                                               + pow(Site::originTransformVector()(k), 2), rPower/2);
            }
        }
    }

    m_potential *= scale;

}

double DiffusionReaction::getSaddleEnergy()
{

    double xs = ((x() + xD())%NX)/2.0;
    double ys = ((y() + yD())%NY)/2.0;
    double zs = ((z() + zD())%NZ)/2.0;

    vector<const Site*> neighborSet;

    for (const Site* site : m_reactionSite->allNeighbors())
    {
        for (const Site* dSite : destination->allNeighbors())
        {
            if (site == dSite)
            {
                neighborSet.push_back(site);
            }
        }
    }

#ifndef NDEBUG

    uint a = neighborSet.size();
    bool b = a == 98 || a == 62 || a == 78;

    if (!b)
    {
        set<Site*> nSet;
        vector<Site*> lSet;
        cout << "THIS SHOULD NEVER OCCUR" << endl;
        cout << destination->allNeighbors().size() << endl;
        cout << m_reactionSite->allNeighbors().size() << endl;

        uint C = 0;
        uint K = 0;
        for (Site* site : m_reactionSite->allNeighbors())
        {
            uint L = 0;
            for (Site* siteD : destination->allNeighbors())
            {
                if (site == siteD)
                {
                    cout << "true " << site << " " << siteD << endl;
                    C++;
                    nSet.insert(site);
                    lSet.push_back(site);
                }

                if (!(siteD == *(destination->allNeighbors().begin() + L)))
                {
                   cout << "fail" << endl;
                   exit(1);
                }
                L++;
            }

            if (!(site == *(m_reactionSite->allNeighbors().begin() + K)))
            {
                cout << "fail" << endl;
                exit(1);
            }

            K++;
        }
        std::set_intersection(m_reactionSite->allNeighbors().begin(), m_reactionSite->allNeighbors().end(),
                              destination->allNeighbors().begin(), destination->allNeighbors().end(),
                              std::inserter(neighborSet, neighborSet.end()));

        cout << C << endl;


        dumpInfo();
        cout << lSet.size() << " " << nSet.size() << endl;
        exit(1);
    }

#endif

    double Esp = 0;

    for (const Site* targetSite : neighborSet)
    {

        if (targetSite->isActive())
        {

            double dx = fabs(xs - targetSite->x());
            double dy = fabs(ys - targetSite->y());
            double dz = fabs(zs - targetSite->z());

            if (dx > Site::nNeighborsLimit())
            {
                dx = NX - dx;
            }

            if (dy > Site::nNeighborsLimit())
            {
                dy = NY - dy;
            }

            if (dz > Site::nNeighborsLimit())
            {
                dz = NZ - dz;
            }


            double r = sqrt(dx*dx + dy*dy + dz*dz);

            assert(r >= 1/2. && "Saddle point is atleast this distance from another site.");
            Esp += scale/pow(r, rPower);
        }
    }

    lastUsedEsp = Esp;

    return Esp;

}

void DiffusionReaction::calcRate()
{

    lastUsedE = m_reactionSite->energy();
    m_rate = m_linearRateScale*exp(-beta*(m_reactionSite->energy()-getSaddleEnergy()));

}

bool DiffusionReaction::isNotBlocked()
{

    return !destination->isActive() && (destination->isSurface() || (destination->nNeighbors() == 1));

}

void DiffusionReaction::execute()
{

    m_reactionSite->deactivate();
    destination->activate();

}

void DiffusionReaction::dumpInfo(int xr, int yr, int zr)
{

    (void) xr;
    (void) yr;
    (void) zr;

    int X, Y, Z;
    m_reactionSite->distanceTo(destination, X, Y, Z);

    assert((x() + NX + X)%NX == destination->x());
    assert((y() + NY + Y)%NY == destination->y());
    assert((z() + NZ + Z)%NZ == destination->z());

    Reaction::dumpInfo(X, Y, Z);

    cout << "Path: " << X << " " << Y << " " << Z << endl;
    cout << "Reaction initiates diffusion to " << endl;
    cout << "{\n";
    destination->dumpInfo(-X, -Y, -Z);
    cout << "\n}" << endl;

}

bool DiffusionReaction::allowedAtSite()
{

    //Diffusion reactions may occur to surfaces
    if (destination->isSurface())
    {
        return true;
    }

    //if were not on a surface, we check if te destination is close to other particles.
    else
    {
        return destination->nNeighbors() == 0;
    }

    return true;

}


double DiffusionReaction::rPower = 0;
double DiffusionReaction::scale = 0;

cube   DiffusionReaction::m_potential;
