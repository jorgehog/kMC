#include "diffusionreaction.h"
#include "../../kmcsolver.h"

#include "../../debugger/debugger.h"


using namespace kMC;

DiffusionReaction::DiffusionReaction(Site * currentSite, Site *destinationSite) :
    Reaction(currentSite, "DiffusionReaction"),
    m_lastUsedEsp(UNSET_ENERGY),
    m_destinationSite(destinationSite)
{

    reactionSite()->distanceTo(destinationSite, path(0), path(1), path(2));

    neighborSetIntersectionPoints = getSaddleOverlapMatrix(path);

    saddleFieldIndices = conv_to<uvec>::from(path + 1);

}

DiffusionReaction::~DiffusionReaction()
{
    m_destinationSite = NULL;
}


void DiffusionReaction::loadConfig(const Setting &setting)
{


    separation = getSurfaceSetting<uint>(setting, "separation");

    allowanceFunction = [] (const DiffusionReaction* reaction)
    {
        return (reaction->destinationSite()->isSurface() ||
                reaction->destinationSite()->nNeighbors() == (reaction->reactionSite()->isActive() ? 1 : 0));
    };




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

    umat overlapBox;
    ivec _path;
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

                _path = {x, y, z};
                overlapBox = getSaddleOverlapMatrix(_path);


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

                            r = sqrt(dx*dx + dy*dy + dz*dz);


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

    if (!Debugger::enabled) return "";

    int X, Y, Z;
    X = 0;
    Y = 0;
    Z = 0;

    stringstream s;

    s << Reaction::getFinalizingDebugMessage();

    const Reaction * lastReaction = Debugger::lastCurrentReaction;

    if (lastReaction != NULL)
    {
        const Site* dest = static_cast<const DiffusionReaction*>(lastReaction)->destinationSite();
        reactionSite()->distanceTo(dest, X, Y, Z);
    }

    s << "\nDestination of last active reaction site marked on current site:\n\n";
    s << reactionSite()->info(X, Y, Z);

    return s.str();
#else
    return "";
#endif
}


void DiffusionReaction::setDirectUpdateFlags(const Site *changedSite)
{

    uint d_maxDistance, r_maxDistance;


    if (rate() == UNSET_RATE || changedSite == reactionSite())
    {
        addUpdateFlag(defaultUpdateFlag);
    }

    else
    {

        r_maxDistance = reactionSite()->maxDistanceTo(changedSite);

        KMCDebugger_Assert(r_maxDistance, !=, 0, "This should be handled by other test.", getFinalizingDebugMessage());

        if (r_maxDistance == 1)
        {
            addUpdateFlag(defaultUpdateFlag);
        }

        else
        {
            d_maxDistance = m_destinationSite->maxDistanceTo(changedSite);

            //if the destination is outsite the interaction cutoff, we can keep the old saddle energy.
            if (d_maxDistance > Site::nNeighborsLimit())
            {
                KMCDebugger_Assert(Site::nNeighborsLimit() + 1, ==,  d_maxDistance);
                addUpdateFlag(updateKeepSaddle);
            }

            else
            {
                addUpdateFlag(defaultUpdateFlag);
            }
        }


        KMCDebugger_AssertBool(!updateFlags().empty(), "Updateflag should not be empty!", info());
    }


}

double DiffusionReaction::getSaddleEnergy()
{

    if (reactionSite()->nNeighborsSum() == 0 || m_destinationSite->nNeighborsSum() == 1)
    {
        return 0;
    }

    Site * targetSite;

    double Esp = 0;

    const cube & saddlePot = m_saddlePotential(saddleFieldIndices(0),
                                               saddleFieldIndices(1),
                                               saddleFieldIndices(2));

    for (uint xn = neighborSetIntersectionPoints(0, 0); xn < neighborSetIntersectionPoints(0, 1); ++xn)
    {
        for (uint yn = neighborSetIntersectionPoints(1, 0); yn < neighborSetIntersectionPoints(1, 1); ++yn)
        {
            for (uint zn = neighborSetIntersectionPoints(2, 0); zn < neighborSetIntersectionPoints(2, 1); ++zn)
            {
                targetSite = reactionSite()->neighborHood(xn, yn, zn);

                if (targetSite == NULL)
                {
                    continue;
                }

                else if (!targetSite->isActive())
                {
                    continue;
                }

                else if (targetSite == reactionSite())
                {
                    continue;
                }

                Esp += saddlePot(xn - neighborSetIntersectionPoints(0, 0),
                                 yn - neighborSetIntersectionPoints(1, 0),
                                 zn - neighborSetIntersectionPoints(2, 0));

                KMCDebugger_Assert(saddlePot(xn - neighborSetIntersectionPoints(0, 0),
                                             yn - neighborSetIntersectionPoints(1, 0),
                                             zn - neighborSetIntersectionPoints(2, 0)),
                                   ==,
                                   getSaddleEnergyContributionFrom(targetSite),
                                   "Mismatch in saddle energy contribution.",
                                   getFinalizingDebugMessage());
            }
        }
    }


    return Esp;

}

double DiffusionReaction::getSaddleEnergyContributionFrom(const Site *site)
{
    int X, Y, Z;

    reactionSite()->distanceTo(site, X, Y, Z);

    return getSaddleEnergyContributionFromNeighborAt(X + Site::nNeighborsLimit(),
                                                     Y + Site::nNeighborsLimit(),
                                                     Z + Site::nNeighborsLimit());

}

double DiffusionReaction::getSaddleEnergyContributionFromNeighborAt(const uint &i, const uint &j, const uint &k)
{
    return m_saddlePotential(saddleFieldIndices(0),
                             saddleFieldIndices(1),
                             saddleFieldIndices(2))
            (i - neighborSetIntersectionPoints(0, 0),
             j - neighborSetIntersectionPoints(1, 0),
             k - neighborSetIntersectionPoints(2, 0));

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

    double newRate = 0;

    if (updateFlag() == defaultUpdateFlag)
    {

        double Esp = getSaddleEnergy();

        newRate = linearRateScale()*exp(-beta()*(reactionSite()->energy()- Esp));

        m_lastUsedEsp = Esp;
    }

    else //m_udateFlag = updateKeepSaddle
    {

        KMCDebugger_Assert(rate(), !=, UNSET_RATE ,"Saddle can't update when the rate has not been calculated.", getFinalizingDebugMessage());
        KMCDebugger_Assert(updateFlag(), ==, updateKeepSaddle, "Errorous updateFlag.", getFinalizingDebugMessage());
        KMCDebugger_Assert(lastUsedEnergy(), !=, UNSET_ENERGY, "energy never calculated before.", getFinalizingDebugMessage());

        newRate = rate()*exp(-beta()*(reactionSite()->energy() - lastUsedEnergy()));

        KMCDebugger_AssertClose(getSaddleEnergy(), m_lastUsedEsp, 1E-10, "Saddle energy was not conserved as assumed by flag. ", getFinalizingDebugMessage());

    }

    setRate(newRate);

}

void DiffusionReaction::execute()
{
    reactionSite()->deactivate();
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

    stringstream s;
    s << Reaction::info(X, Y, Z, "D");

    s << "Reaction initiates diffusion to\n\n";

    s << m_destinationSite->info(-X, -Y, -Z, "O");

    s << "\nPath: " << X << " " << Y << " " << Z << endl;

    string full_string = s.str();

    return full_string;

}

bool DiffusionReaction::isAllowed() const
{

    return !m_destinationSite->isActive() && allowanceFunction(this);

}


double        DiffusionReaction::rPower;
double        DiffusionReaction::scale;

uint          DiffusionReaction::separation;

cube          DiffusionReaction::m_potential;
field<cube>   DiffusionReaction::m_saddlePotential;

function<bool (const DiffusionReaction*)> DiffusionReaction::allowanceFunction;
