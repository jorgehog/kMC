#include "diffusionreaction.h"
#include "../../kmcsolver.h"

DiffusionReaction::DiffusionReaction(Site *destination) :
    Reaction("DiffusionReaction"),
    destination(destination),
    updateFlag(updateFull),
    energyShift(0)
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

    //rescale the potential to avoid exploding rates for some choices of parameters.
    //    scale = 1.0/accu(m_potential);

    m_potential *= scale;

}


void DiffusionReaction::setUpdateFlags(const Site *changedSite, uint i, uint j, uint k, uint level, double dE)
{

    (void) i;
    (void) j;
    (void) k;

    if (level == 0 || m_rate == UNSET_RATE)
    {
        updateFlag = updateFull;
    }

    else if (destination->maxDistanceTo(changedSite) == Site::nNeighborsLimit() + 1)
    {
        energyShift = dE;
        updateFlag = updateNoSaddle;
    }

    else
    {
        assert(level == Site::nNeighborsLimit() - 1);
        updateFlag = updateFull;
    }

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
            if (site == dSite && site->isActive())
            {
                neighborSet.push_back(site);
            }
        }
    }

    double Esp = 0;

    for (const Site* targetSite : neighborSet)
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

    if (lastUsedEsp == Esp)
    {
        counter++;
    }
    total++;

    lastUsedEsp = Esp;

    return Esp;

}

void DiffusionReaction::calcRate()
{

    if (updateFlag == updateNoSaddle)
    {
        assert(m_rate != UNSET_RATE);
        assert(energyShift != 0);

        m_rate *= exp(-beta*energyShift);

#ifndef NDEBUG

        if (fabs(lastUsedE - reactionSite()->energy() + energyShift) > 1E-15)
        {
            cout << "ENERGY SHIFT MISMATCH" << endl;
            cout << lastUsedE - reactionSite()->energy() <<endl;
            cout << -energyShift << endl;
            cout << lastUsedE - reactionSite()->energy() + energyShift << endl;
            exit(1);
        }
        double newRate = m_rate;
        updateFlag = updateFull;
        calcRate();
        if (fabs(newRate - m_rate) > 1E-15)
        {
            cout << "SOMETHING WENT WRONG IN UPDATEING ALG" << endl;
            cout << newRate << "  " << m_rate << endl;
            cout << newRate - m_rate << endl;
            exit(1);
        }

#endif

    }

    else
    {
        m_rate = m_linearRateScale*exp(-beta*(m_reactionSite->energy()-getSaddleEnergy()));
    }

    lastUsedE = m_reactionSite->energy();

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

uint   DiffusionReaction::total = 0;
uint   DiffusionReaction::counter = 0;
