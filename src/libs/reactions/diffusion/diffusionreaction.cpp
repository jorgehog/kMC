#include "diffusionreaction.h"
#include "../../kmcsolver.h"

#include "../../debugger/kmcdebugger.h"

DiffusionReaction::DiffusionReaction(Site *destination) :
    Reaction("DiffusionReaction"),
    lastUsedEnergy(0),
    lastUsedEsp(0),
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

    //rescale the potential to avoid exploding rates for some choices of parameters.
    //    scale = 1.0/accu(m_potential);

    m_potential *= scale;

}
/*
 * 31% 64% old | 5% 47 % new
 */

void DiffusionReaction::setUpdateFlags(const Site *changedSite, uint level)
{

    if (m_rate == UNSET_RATE)
    {
        m_updateFlags.insert(defaultUpdateFlag);
    }

    else if (changedSite->isSurface())
    {
        if (!destination->isSurface())
        {
            m_updateFlags.insert(updateKeepSaddle);
        }
    }

    else if (level == 0)
    {
        m_updateFlags.insert(defaultUpdateFlag);
    }

    //if the destination is outsite the interaction cutoff, we can keep the old saddle energy.
    else if (destination->maxDistanceTo(changedSite) == Site::nNeighborsLimit() + 1)
    {
        m_updateFlags.insert(updateKeepSaddle);
    }

    else
    {
        assert(level == Site::nNeighborsLimit() - 1);
        m_updateFlags.insert(defaultUpdateFlag);
    }

}

double DiffusionReaction::getSaddleEnergy()
{

    timer.tic();

    double xs = ((x() + xD())%NX)/2.0;
    double ys = ((y() + yD())%NY)/2.0;
    double zs = ((z() + zD())%NZ)/2.0;

    vector<const Site*> neighborSet;
    bool empty = true;

    for (const Site* site : m_reactionSite->allNeighbors())
    {
        for (const Site* dSite : destination->allNeighbors())
        {
            if (site == dSite && site->isActive())
            {
                neighborSet.push_back(site);
                empty = false;
            }
        }
    }

    if (empty)
    {
        return 0;
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

    bool sameSetup = true;

    if (lastSetup.size() != neighborSet.size())
    {
        sameSetup = false;
    }
    else
    {
        for (const Site* s: neighborSet)
        {
            bool isIn = false;
            for (const Site* slast : lastSetup)
            {
                if (s == slast)
                {
                    isIn = true;
                    break;
                }
            }
            if (!isIn)
            {
                sameSetup = false;
                break;
            }
        }
    }

    if (fabs(lastUsedEsp - Esp) < 1E-10)
    {
        if (sameSetup)
        {

            cout << "exactly same setup calculated saddle twice..should be flagged" << endl;

            KMCDebugger_DumpFullTrace(true, info());

            exit(1);
        }
        counterEqSP++;

    }

    totalSP++;

    lastUsedEsp = Esp;

    for (const Site * s: neighborSet)
    {
        lastSetup.insert(s);
    }

    totalTime += timer.toc();

    return Esp;

}

void DiffusionReaction::calcRate()
{

    counterAllRate++;

    switch (m_updateFlag) {
    case noUpdate:

        assert(m_reactionSite->energy() == lastUsedEnergy);

        return;

    case defaultUpdateFlag:

        m_rate = m_linearRateScale*exp(-beta*(m_reactionSite->energy()-getSaddleEnergy()));

        break;

    case updateKeepSaddle:

        m_rate *= exp(-beta*(reactionSite()->energy() - lastUsedEnergy));

        break;

    default:

        cout << "Unknown update flag: " << m_updateFlag << endl;
        exit(1);

        break;
    }

    lastUsedEnergy = m_reactionSite->energy();

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

const string DiffusionReaction::info(int xr, int yr, int zr, string desc) const
{

    (void) xr;
    (void) yr;
    (void) zr;
    (void) desc;

    int X, Y, Z;
    m_reactionSite->distanceTo(destination, X, Y, Z);

    assert((x() + NX + X)%NX == destination->x());
    assert((y() + NY + Y)%NY == destination->y());
    assert((z() + NZ + Z)%NZ == destination->z());

    stringstream s;
    s << Reaction::info(X, Y, Z, "D");

    s << "Path: " << X << " " << Y << " " << Z << endl;
    s << "Reaction initiates diffusion to " << endl;
    s << "{\n";
    s << destination->info(-X, -Y, -Z, "O");
    s << "\n}" << endl;

    return s.str();

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


double DiffusionReaction::rPower ;
double DiffusionReaction::scale;

cube   DiffusionReaction::m_potential;

uint   DiffusionReaction::totalSP = 0;
uint   DiffusionReaction::counterEqSP = 0;
uint   DiffusionReaction::counterAllRate = 0;
double DiffusionReaction::totalTime = 0;
wall_clock DiffusionReaction::timer;
