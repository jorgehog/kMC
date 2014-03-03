#include "diffusionreaction.h"
#include "../../kmcsolver.h"

#include "../../debugger/kmcdebugger.h"

DiffusionReaction::DiffusionReaction(Site * currentSite, Site *destinationSite) :
    Reaction("DiffusionReaction"),
    lastUsedEnergy(UNSET_ENERGY),
    lastUsedEsp(UNSET_ENERGY),
    m_destinationSite(destinationSite)
{

    setSite(currentSite);

    m_reactionSite->distanceTo(destinationSite, path(0), path(1), path(2));

    rSaddle = conv_to<vec>::from(path)/2;

    neighborSetIntersectionPoints = getSaddleOverlapMatrix(path);

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

    m_saddlePotential.set_size(3, 3, 3);

    uint i, j, k;
    double dx, dy, dz, r;

    for (int x = -1; x <= 1; ++x)
    {
        for (int y = -1; y <= 1; ++y)
        {
            for (int z = -1; z <= 1; ++z)
            {
                if ((x == 0) && (y == 0) && (z == 0))
                {
                    continue;
                }

                i = x + 1;
                j = y + 1;
                k = z + 1;

                umat overlapBox = getSaddleOverlapMatrix({x, y, z});


                m_saddlePotential(i, j, k).set_size(overlapBox(0, 1) - overlapBox(0, 0),
                                                    overlapBox(1, 1) - overlapBox(1, 0),
                                                    overlapBox(2, 1) - overlapBox(2, 0));

                KMCDebugger_Assert(m_saddlePotential.n_elem, !=, 0, "illegal box size.");
                KMCDebugger_Assert(m_saddlePotential.n_elem, !=, arma::prod(overlapBox.col(1)-overlapBox.col(0)), "saddle mat fail.");
//                cout << " --------  " << x << " " << y << " " << z << "  ------ " << endl;
                for (uint xn = overlapBox(0, 0); xn < overlapBox(0, 1); ++xn)
                {

                    dx = x/2.0 + (int)xn - (int)Site::nNeighborsLimit();

                    for (uint yn = overlapBox(1, 0); yn < overlapBox(1, 1); ++yn)
                    {

                        dy = y/2.0 + (int)yn - (int)Site::nNeighborsLimit();

                        for (uint zn = overlapBox(2, 0); zn < overlapBox(2, 1); ++zn)
                        {

                            if (xn == yn && yn == zn && zn == Site::nNeighborsLimit())
                            {
                                continue;
                            }

                            dz = z/2.0 + (int)zn - (int)Site::nNeighborsLimit();

                            r = sqrt(dx*dx + dy*dy + dz*dz);

//                            cout << xn << " " << yn << " " << zn << " " << dx << " " << dy << " " << dz << endl;

                            m_saddlePotential(i, j, k)(xn - overlapBox(0, 0),
                                                       yn - overlapBox(1, 0),
                                                       zn - overlapBox(2, 0)) = 1.0/pow(r, rPower);

                        }
                    }
                }

                m_saddlePotential(i, j, k) *= scale;

            }
        }
    }




    //rescale the potential to avoid exploding rates for some choices of parameters.
    //    scale = 1.0/accu(m_potential);
    m_potential *= scale;

}

const uint &DiffusionReaction::xD() const
{
    return m_destinationSite->x();
}

const uint &DiffusionReaction::yD() const
{
    return m_destinationSite->y();
}

const uint &DiffusionReaction::zD() const
{
    return m_destinationSite->z();
}

string DiffusionReaction::getFinalizingDebugMessage() const
{
#ifndef KMC_NO_DEBUG

    if (!KMCDebugger::enabled) return "";

    int X, Y, Z;
    X = 0;
    Y = 0;
    Z = 0;

    stringstream s;

    s << Reaction::getFinalizingDebugMessage();

    const Reaction * lastReaction = KMCDebugger::lastCurrentReaction;

    if (lastReaction != NULL)
    {
        const Site* dest = static_cast<const DiffusionReaction*>(lastReaction)->destinationSite();
        m_reactionSite->distanceTo(dest, X, Y, Z);
    }

    s << "\nDestination of last active reaction site marked on current site:\n\n";
    s << m_reactionSite->info(X, Y, Z);

    return s.str();
#else
    return "";
#endif
}


void DiffusionReaction::setDirectUpdateFlags(const Site *changedSite)
{

    uint d_maxDistance, r_maxDistance;


    if (m_rate == UNSET_RATE || changedSite == m_reactionSite)
    {
        setImplicitUpdateFlags();
    }

    else
    {

        r_maxDistance = reactionSite()->maxDistanceTo(changedSite);


        if (r_maxDistance == 1)
        {
            setImplicitUpdateFlags();
        }

        else
        {
            d_maxDistance = m_destinationSite->maxDistanceTo(changedSite);

            //if the destination is outsite the interaction cutoff, we can keep the old saddle energy.
            if (d_maxDistance > Site::nNeighborsLimit())
            {
                KMCDebugger_Assert(Site::nNeighborsLimit() + 1, ==,  d_maxDistance);
                m_updateFlags.insert(updateKeepSaddle);
            }

            else
            {
                KMCDebugger_Assert(r_maxDistance, ==, Site::nNeighborsLimit());
                setImplicitUpdateFlags();
            }
        }


        KMCDebugger_AssertBool(!m_updateFlags.empty(), "Updateflag should not be empty!", info());
    }


}

double DiffusionReaction::getSaddleEnergy()
{

    if (m_reactionSite->nNeighborsSum() == 0 || m_destinationSite->nNeighborsSum() == 1)
    {
        return 0;
    }


    timer.tic();


    Site * targetSite;

    double xs, ys, zs;

    xs = x() + 0.5*path(0);
    ys = y() + 0.5*path(1);
    zs = z() + 0.5*path(2);

    double Esp = 0;
    double Esp2 = 0;

    cout << path.t();
    for (uint xn = neighborSetIntersectionPoints(0, 0); xn < neighborSetIntersectionPoints(0, 1); ++xn)
    {
        for (uint yn = neighborSetIntersectionPoints(1, 0); yn < neighborSetIntersectionPoints(1, 1); ++yn)
        {
            for (uint zn = neighborSetIntersectionPoints(2, 0); zn < neighborSetIntersectionPoints(2, 1); ++zn)
            {
                targetSite = m_reactionSite->neighborHood(xn, yn, zn);

                if (!targetSite->isActive())
                {
                    continue;
                }

                else if (targetSite == m_reactionSite)
                {
                    continue;
                }

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

                cout << xn << " " << yn << " " << zn << " " << dx << " " << dy << " " << dz << endl;

                assert(r >= 1/2. && "Saddle point is atleast this distance from another site.");
                Esp += scale/pow(r, rPower);


                const cube & saddlePot = m_saddlePotential(path(0) + 1, path(1) + 1, path(2) + 1);

                Esp2 += saddlePot(xn - neighborSetIntersectionPoints(0, 0),
                                  yn - neighborSetIntersectionPoints(1, 0),
                                  zn - neighborSetIntersectionPoints(2, 0));

                KMCDebugger_Assert(Esp, ==, Esp2);
//                if (Esp == Esp2)
//                {
//                    cout << "WIN"<< endl;
//                } else {
//                    cout << "FAIL" << endl;
//                }

            }
        }
    }

    cout << "######################" << endl;

    totalTime += timer.toc();

    return Esp;

}

umat::fixed<3, 2> DiffusionReaction::getSaddleOverlapMatrix(const ivec & relCoor)
{

    umat::fixed<3, 2> overlap;

    for (uint xyz = 0; xyz < 3; ++xyz)
    {
        if (relCoor(xyz) == 1)
        {

            overlap(xyz, 0) = 1;
            overlap(xyz, 1) = Site::neighborhoodLength();

        }

        else if (relCoor(xyz) == -1)
        {

            overlap(xyz, 0) = 0;
            overlap(xyz, 1) = Site::neighborhoodLength() - 1;
        }

        else
        {

            overlap(xyz, 0) = 0;
            overlap(xyz, 1) = Site::neighborhoodLength();

            KMCDebugger_Assert(overlap(xyz), ==, 0, "There should be no other option.");
        }
    }

    return overlap;
}

void DiffusionReaction::calcRate()
{

    counterAllRate++;

    if (m_updateFlag == defaultUpdateFlag)
    {

        double Esp = getSaddleEnergy();

        m_rate = m_linearRateScale*exp(-beta*(m_reactionSite->energy()- Esp));

        if (Esp != 0)
        {

            if (fabs(lastUsedEsp - Esp) < 1E-10)
            {
                //                KMCDebugger_AssertBool(false, "should never happen.", getFinalizingDebugMessage());
                counterEqSP++;
            }

        }

        totalSP++;
        lastUsedEsp = Esp;
    }

    else //m_udateFlag = updateKeepSaddle
    {

        KMCDebugger_Assert(m_rate, !=, UNSET_RATE ,"Saddle can't update when the rate has not been calculated.", getFinalizingDebugMessage());
        KMCDebugger_Assert(m_updateFlag, ==, updateKeepSaddle, "Errorous updateFlag.", getFinalizingDebugMessage());
        KMCDebugger_Assert(lastUsedEnergy, !=, UNSET_ENERGY, "energy never calculated before.", getFinalizingDebugMessage());

        m_rate *= exp(-beta*(reactionSite()->energy() - lastUsedEnergy));

        KMCDebugger_AssertClose(getSaddleEnergy(), lastUsedEsp, 1E-10, "Saddle energy was not conserved as assumed by flag. ", getFinalizingDebugMessage());

    }

    lastUsedEnergy = m_reactionSite->energy();

}

bool DiffusionReaction::isNotBlocked() const
{

    return !m_destinationSite->isActive() && (m_destinationSite->isSurface() || (m_destinationSite->nNeighbors() == 1));

}

void DiffusionReaction::execute()
{
    m_reactionSite->deactivate();
    m_destinationSite->activate();

}

const string DiffusionReaction::info(int xr, int yr, int zr, string desc) const
{

    (void) xr;
    (void) yr;
    (void) zr;
    (void) desc;

    int X = path(0);
    int Y = path(1);
    int Z = path(2);

    assert((x() + NX + X)%NX == m_destinationSite->x());
    assert((y() + NY + Y)%NY == m_destinationSite->y());
    assert((z() + NZ + Z)%NZ == m_destinationSite->z());

    stringstream s;
    s << Reaction::info(X, Y, Z, "D");

    s << "Reaction initiates diffusion to\n\n";

    s << m_destinationSite->info(-X, -Y, -Z, "O");

    s << "\nPath: " << X << " " << Y << " " << Z << endl;

    return s.str();

}

bool DiffusionReaction::allowedAtSite()
{

    //Diffusion reactions may occur to surfaces
    if (m_destinationSite->isSurface())
    {
        return true;
    }

    //if were not on a surface, we check if te destination is close to other particles.
    else
    {
        return m_destinationSite->nNeighbors() == 0;
    }

    return true;

}

const
double        DiffusionReaction::UNSET_ENERGY = -1;

double        DiffusionReaction::rPower;
double        DiffusionReaction::scale;

cube          DiffusionReaction::m_potential;
field<cube>   DiffusionReaction::m_saddlePotential;

uint          DiffusionReaction::totalSP = 0;
uint          DiffusionReaction::counterEqSP = 0;
uint          DiffusionReaction::counterAllRate = 0;
double        DiffusionReaction::totalTime = 0;
wall_clock    DiffusionReaction::timer;
