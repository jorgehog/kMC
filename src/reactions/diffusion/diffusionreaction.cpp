#include "diffusionreaction.h"

#include "../../kmcsolver.h"

#include "../../boundary/boundary.h"

#include "../../soluteparticle.h"

#include "../../debugger/debugger.h"


using namespace kMC;

DiffusionReaction::DiffusionReaction(SoluteParticle *reactant, int dx, int dy, int dz) :
    Reaction(reactant),
    m_lastUsedEsp(UNSET_ENERGY)
{

    m_path[0] = dx;
    m_path[1] = dy;
    m_path[2] = dz;

    saddleFieldIndices[0] = dx + 1;
    saddleFieldIndices[1] = dy + 1;
    saddleFieldIndices[2] = dz + 1;

}

DiffusionReaction::~DiffusionReaction()
{

}


void DiffusionReaction::loadConfig(const Setting &setting)
{

    m_rPower = getSurfaceSetting<double>(setting, "rPower");
    m_scale  = getSurfaceSetting<double>(setting, "scale");

}

uint DiffusionReaction::xD() const
{
    Boundary::setupCurrentBoundary(x(), 0);
    return Boundary::currentBoundaries(0)->transformCoordinate(x() + m_path[0]);
}

uint DiffusionReaction::yD() const
{
    Boundary::setupCurrentBoundary(y(), 1);
    return Boundary::currentBoundaries(1)->transformCoordinate(y() + m_path[1]);
}

uint DiffusionReaction::zD() const
{
    Boundary::setupCurrentBoundary(z(), 2);
    return Boundary::currentBoundaries(2)->transformCoordinate(z() + m_path[2]);
}

string DiffusionReaction::getFinalizingDebugMessage() const
{
#ifndef KMC_NO_DEBUG

    if (!Debugger::enabled) return "";

    int X, Y, Z;
    X = 0;
    Y = 0;
    Z = 0;

    stringstream s;

    s << Reaction::getFinalizingDebugMessage();


    if (!Debugger::lastCurrentReaction->isType("DiffusionReaction"))
    {
        return s.str();
    }

    const DiffusionReaction *lastReaction = (DiffusionReaction*)Debugger::lastCurrentReaction;

    if (lastReaction != NULL)
    {
        X = lastReaction->path(0);
        Y = lastReaction->path(1);
        Z = lastReaction->path(2);
    }

    s << "\nDestination of last active reaction site marked on current site:\n\n";
    s << reactant()->info(X, Y, Z);

    return s.str();
#else
    return "";
#endif
}

void DiffusionReaction::setupPotential()
{

    KMCDebugger_Assert(m_scale, !=, 0, "Potential parameters not set.");

    m_potential.reset();
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

                m_potential(i, j, k) = 1.0/std::pow(Site::originTransformVector(i)*Site::originTransformVector(i)
                                                    + Site::originTransformVector(j)*Site::originTransformVector(j)
                                                    + Site::originTransformVector(k)*Site::originTransformVector(k)
                                                    , m_rPower/2);
            }
        }
    }

    m_saddlePotential.reset_objects();
    m_saddlePotential.reset();
    m_saddlePotential.set_size(3, 3, 3);

    neighborSetIntersectionPoints.reset_objects();
    neighborSetIntersectionPoints.reset();
    neighborSetIntersectionPoints.set_size(3, 3, 3);

    umat::fixed<3, 2> overlapBox;
    ivec _path;
    uint i, j, k;
    double dx, dy, dz, r2;

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

                _path = {x, y, z};
                overlapBox = makeSaddleOverlapMatrix(_path);

                neighborSetIntersectionPoints(i, j, k) = overlapBox;

                m_saddlePotential(i, j, k).set_size(overlapBox(0, 1) - overlapBox(0, 0),
                                                    overlapBox(1, 1) - overlapBox(1, 0),
                                                    overlapBox(2, 1) - overlapBox(2, 0));

                KMCDebugger_Assert(m_saddlePotential.n_elem, !=, 0, "illegal box size.");
                KMCDebugger_Assert(m_saddlePotential.n_elem, !=, arma::prod(overlapBox.col(1)-overlapBox.col(0)), "saddle mat fail.");

                for (uint xn = overlapBox(0, 0); xn < overlapBox(0, 1); ++xn)
                {

                    dx = x/2.0 - (int)xn + (int)Site::nNeighborsLimit();

                    for (uint yn = overlapBox(1, 0); yn < overlapBox(1, 1); ++yn)
                    {

                        dy = y/2.0 - (int)yn + (int)Site::nNeighborsLimit();

                        for (uint zn = overlapBox(2, 0); zn < overlapBox(2, 1); ++zn)
                        {

                            if (xn == yn && yn == zn && zn == Site::nNeighborsLimit())
                            {
                                continue;
                            }

                            dz = z/2.0 - (int)zn + (int)Site::nNeighborsLimit();

                            r2 = dx*dx + dy*dy + dz*dz;


                            m_saddlePotential(i, j, k)(xn - overlapBox(0, 0),
                                                       yn - overlapBox(1, 0),
                                                       zn - overlapBox(2, 0)) = 1.0/pow(r2, m_rPower/2);

                        }
                    }
                }

                m_saddlePotential(i, j, k) *= m_scale;

            }
        }
    }


    m_potential *= m_scale;

}

Site *DiffusionReaction::destinationSite() const
{
    return Site::neighborhood(x(), y(), z(), m_path[0], m_path[1], m_path[2]);
}


void DiffusionReaction::setDirectUpdateFlags(const SoluteParticle *changedReactant, const uint level)
{

    KMCDebugger_Assert(changedReactant, !=, reactant());


    if (rate() == UNSET_RATE)
    {
        registerUpdateFlag(defaultUpdateFlag);
    }

    else
    {

        //        r_maxDistance = site()->maxDistanceTo(changedReactant->site());

        if (level == 0)
        {
            registerUpdateFlag(defaultUpdateFlag);
        }

        else
        {

            uint d_maxDistance = reactant()->maxDistanceTo(changedReactant);

            //if the destination is outsite the interaction cutoff, we can keep the old saddle energy.
            if (d_maxDistance > Site::nNeighborsLimit())
            {
                KMCDebugger_Assert(Site::nNeighborsLimit() + 1, ==,  d_maxDistance);
                registerUpdateFlag(updateKeepSaddle);
            }

            else
            {
                registerUpdateFlag(defaultUpdateFlag);
            }
        }

    }


}

double DiffusionReaction::getSaddleEnergy()
{

    if (reactant()->nNeighborsSum() == 0)
    {
        return 0;
    }

    Site * targetSite;

    double Esp = 0;

    const cube & saddlePot = m_saddlePotential(saddleFieldIndices[0],
            saddleFieldIndices[1],
            saddleFieldIndices[2]);

    const umat::fixed<3, 2> & myIntersectionPoints = neighborSetIntersectionPoints(saddleFieldIndices[0],
            saddleFieldIndices[1],
            saddleFieldIndices[2]);

    for (uint xn = myIntersectionPoints(0, 0); xn < myIntersectionPoints(0, 1); ++xn)
    {
        for (uint yn = myIntersectionPoints(1, 0); yn < myIntersectionPoints(1, 1); ++yn)
        {
            for (uint zn = myIntersectionPoints(2, 0); zn < myIntersectionPoints(2, 1); ++zn)
            {

                if (xn == yn && yn == zn && zn == Site::nNeighborsLimit())
                {
                    continue;
                }

                targetSite = Site::neighborhood_fromIndex(x(), y(), z(), xn, yn, zn);

                if (targetSite == NULL)
                {
                    continue;
                }

                else if (!targetSite->isActive())
                {
                    continue;
                }

                Esp += saddlePot(xn - myIntersectionPoints(0, 0),
                                 yn - myIntersectionPoints(1, 0),
                                 zn - myIntersectionPoints(2, 0));

                KMCDebugger_Assert(saddlePot(xn - myIntersectionPoints(0, 0),
                                             yn - myIntersectionPoints(1, 0),
                                             zn - myIntersectionPoints(2, 0)),
                                   ==,
                                   getSaddleEnergyContributionFrom(targetSite->associatedParticle()),
                                   "Mismatch in saddle energy contribution.",
                                   getFinalizingDebugMessage());
            }
        }
    }


    return Esp;

}

double DiffusionReaction::getSaddleEnergyContributionFrom(const SoluteParticle *particle)
{
    int X, Y, Z;

    reactant()->distanceTo(particle, X, Y, Z);

    return getSaddleEnergyContributionFromNeighborAt(X + Site::nNeighborsLimit(),
                                                     Y + Site::nNeighborsLimit(),
                                                     Z + Site::nNeighborsLimit());

}

double DiffusionReaction::getSaddleEnergyContributionFromNeighborAt(const uint &i, const uint &j, const uint &k)
{
    return m_saddlePotential(saddleFieldIndices[0], saddleFieldIndices[1], saddleFieldIndices[2])
            (i - neighborSetIntersectionPoints(saddleFieldIndices[0], saddleFieldIndices[1], saddleFieldIndices[2])(0, 0),
            j - neighborSetIntersectionPoints(saddleFieldIndices[0], saddleFieldIndices[1], saddleFieldIndices[2])(1, 0),
            k - neighborSetIntersectionPoints(saddleFieldIndices[0], saddleFieldIndices[1], saddleFieldIndices[2])(2, 0));
}


umat::fixed<3, 2> DiffusionReaction::makeSaddleOverlapMatrix(const ivec & relCoor)
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

    double newRate = 0;

    KMCDebugger_Assert(updateFlag(), !=, UNSET_UPDATE_FLAG);

    if (updateFlag() == defaultUpdateFlag)
    {

        double Esp = getSaddleEnergy();

        newRate = linearRateScale()*std::exp(-beta()*(reactant()->energy()- Esp));

        m_lastUsedEsp = Esp;
    }

    else //m_udateFlag = updateKeepSaddle
    {

        KMCDebugger_Assert(rate(), !=, UNSET_RATE ,"Saddle can't update when the rate has not been calculated.", getFinalizingDebugMessage());
        KMCDebugger_Assert(updateFlag(), ==, updateKeepSaddle, "Errorous updateFlag.", getFinalizingDebugMessage());
        KMCDebugger_Assert(lastUsedEnergy(), !=, UNSET_ENERGY, "energy never calculated before.", getFinalizingDebugMessage());

        newRate = rate()*std::exp(-beta()*(reactant()->energy() - lastUsedEnergy()));

        KMCDebugger_AssertClose(getSaddleEnergy(), m_lastUsedEsp, 1E-10, "Saddle energy was not conserved as assumed by flag. ", getFinalizingDebugMessage());

    }

    setRate(newRate);

}

void DiffusionReaction::execute()
{
    reactant()->changePosition(Site::boundaries(0, 0)->transformCoordinate(x() + m_path[0]),
                               Site::boundaries(1, 0)->transformCoordinate(y() + m_path[1]),
                               Site::boundaries(2, 0)->transformCoordinate(z() + m_path[2]));
}

const string DiffusionReaction::info(int xr, int yr, int zr, string desc) const
{

    (void) xr;
    (void) yr;
    (void) zr;
    (void) desc;

    stringstream s;
    s << Reaction::info(m_path[0], m_path[1], m_path[2], "D");

    s << "Reaction initiates diffusion to\n\n";

    if (destinationSite() == NULL)
    {
        s << "BOUNDARY";
    }

    else
    {
        s << destinationSite()->info(-m_path[0], -m_path[1], -m_path[2], "O");
    }

    s << "\nPath: " << m_path[0] << " " << m_path[1] << " " << m_path[2] << endl;

    string full_string = s.str();

    return full_string;

}

bool DiffusionReaction::isAllowed() const
{

    if (destinationSite() == NULL)
    {
        return false;
    }

    return !destinationSite()->isActive();
}

void DiffusionReaction::reset()
{

    Reaction::reset();

    m_lastUsedEsp = UNSET_ENERGY;

}

void DiffusionReaction::registerBetaChange(const double newBeta)
{
    (void) newBeta;

    KMCDebugger_AssertBool(rate() != UNSET_RATE, "Beta should not be changed untill rates have been calculated.");

    _setRate(rate() * m_betaChangeScaleFactor);

}


double        DiffusionReaction::m_rPower = 1.0;
double        DiffusionReaction::m_scale  = 1.0;

double        DiffusionReaction::m_betaChangeScaleFactor;

cube          DiffusionReaction::m_potential;
field<cube>   DiffusionReaction::m_saddlePotential;
field<umat::fixed<3, 2> >
DiffusionReaction::neighborSetIntersectionPoints;
