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

    m_saddleFieldIndices[0] = dx + 1;
    m_saddleFieldIndices[1] = dy + 1;
    m_saddleFieldIndices[2] = dz + 1;

}

DiffusionReaction::~DiffusionReaction()
{

}


void DiffusionReaction::loadConfig(const Setting &setting)
{

    const Setting & rPowerSettings = getSetting(setting, "rPowers");
    const Setting & strengthSettings = getSetting(setting, "strengths");

    if (rPowerSettings.getLength() != strengthSettings.getLength())
    {
        KMCSolver::exit("Mismatch in diffusionreaction species parameter lengths.");
    }

    SoluteParticle::nSpecies(rPowerSettings.getLength());

    vector<double> rPowers;
    vector<double> strengths;

    for (uint i = 0; i < SoluteParticle::nSpecies(); ++i)
    {
        rPowers.push_back(rPowerSettings[i]);
        strengths.push_back(strengthSettings[i]);
    }

    setPotentialParameters(rPowers, strengths, false);

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

void DiffusionReaction::setPotentialParameters(const vector<double> &rPowers, const vector<double> &strenghts, bool setup)
{

    if (rPowers.size() != strenghts.size())
    {
        KMCSolver::exit("mismatch in species mixing data.");
    }

    SoluteParticle::nSpecies(rPowers.size());

    m_rPowers.reset();
    m_rPowers.set_size(SoluteParticle::nSpecies(), SoluteParticle::nSpecies());

    m_strengths.reset();
    m_strengths.set_size(SoluteParticle::nSpecies(), SoluteParticle::nSpecies());

    for (uint speciesA = 0; speciesA < SoluteParticle::nSpecies(); ++speciesA)
    {
        for (uint speciesB = 0; speciesB < SoluteParticle::nSpecies(); ++speciesB)
        {
            //geometric mean for powers and arithmetic mean for strengths
            m_rPowers(speciesA, speciesB)   = std::sqrt(rPowers.at(speciesA)*rPowers.at(speciesB));
            m_strengths(speciesA, speciesB) = (strenghts.at(speciesA) + strenghts.at(speciesB))/2;
        }
    }

    if (setup)
    {
        setupPotential();
    }

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

    using std::pow;

    uint i, j, k;
    double dx, dy, dz, r2;
    double rPower, strength;

    umat::fixed<3, 2> overlapBox;

    m_potential.reset_objects();
    m_potential.reset();
    m_potential.set_size(Site::neighborhoodLength(),
                         Site::neighborhoodLength(),
                         Site::neighborhoodLength());

    m_saddlePotential.reset_objects();
    m_saddlePotential.reset();
    m_saddlePotential.set_size(3, 3, 3);

    m_neighborSetIntersectionPoints.reset_objects();
    m_neighborSetIntersectionPoints.reset();
    m_neighborSetIntersectionPoints.set_size(3, 3, 3);


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

                m_neighborSetIntersectionPoints(i, j, k) = makeSaddleOverlapMatrix({x, y, z});
            }
        }
    }


    for (uint i = 0; i < Site::neighborhoodLength(); ++i)
    {
        for (uint j = 0; j < Site::neighborhoodLength(); ++j)
        {
            for (uint k = 0; k < Site::neighborhoodLength(); ++k)
            {

                m_potential(i, j, k).set_size(SoluteParticle::nSpecies(), SoluteParticle::nSpecies());

                for (uint speciesA = 0; speciesA < SoluteParticle::nSpecies(); ++speciesA)
                {
                    for (uint speciesB = 0; speciesB < SoluteParticle::nSpecies(); ++speciesB)
                    {

                        if (i == Site::nNeighborsLimit() && j == Site::nNeighborsLimit() && k == Site::nNeighborsLimit())
                        {
                            m_potential(i, j, k)(speciesA, speciesB) = numeric_limits<double>::max();
                            continue;
                        }
                        rPower = m_rPowers(speciesA, speciesB);
                        strength = m_strengths(speciesA, speciesB);

                        m_potential(i, j, k)(speciesA, speciesB) = strength/std::pow(Site::originTransformVector(i)*Site::originTransformVector(i)
                                                                                   + Site::originTransformVector(j)*Site::originTransformVector(j)
                                                                                   + Site::originTransformVector(k)*Site::originTransformVector(k),
                                                                                     rPower/2);
                    }
                }
            }
        }
    }


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

                overlapBox = m_neighborSetIntersectionPoints(i, j, k);

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


                            m_saddlePotential(i, j, k)(xn - overlapBox(0, 0),
                                                       yn - overlapBox(1, 0),
                                                       zn - overlapBox(2, 0)).set_size(SoluteParticle::nSpecies(), SoluteParticle::nSpecies());

                            dz = z/2.0 - (int)zn + (int)Site::nNeighborsLimit();

                            r2 = dx*dx + dy*dy + dz*dz;

                            for (uint speciesA = 0; speciesA < SoluteParticle::nSpecies(); ++speciesA)
                            {
                                for (uint speciesB = 0; speciesB < SoluteParticle::nSpecies(); ++speciesB)
                                {

                                    rPower = m_rPowers(speciesA, speciesB);
                                    strength = m_strengths(speciesA, speciesB);

                                    m_saddlePotential(i, j, k)(xn - overlapBox(0, 0),
                                                               yn - overlapBox(1, 0),
                                                               zn - overlapBox(2, 0))(speciesA, speciesB) = strength/pow(r2, rPower/2);
                                }
                            }

                        }
                    }
                }

            }
        }
    }

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

    Site *targetSite;
    SoluteParticle *neighbor;

    double Esp = 0;

    const field<mat> & saddlePot = m_saddlePotential(m_saddleFieldIndices[0],
            m_saddleFieldIndices[1],
            m_saddleFieldIndices[2]);

    const umat::fixed<3, 2> & myIntersectionPoints = m_neighborSetIntersectionPoints(m_saddleFieldIndices[0],
            m_saddleFieldIndices[1],
            m_saddleFieldIndices[2]);

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

                neighbor = targetSite->associatedParticle();

                Esp += saddlePot(xn - myIntersectionPoints(0, 0),
                                 yn - myIntersectionPoints(1, 0),
                                 zn - myIntersectionPoints(2, 0))(reactant()->species(), neighbor->species());

                KMCDebugger_Assert(saddlePot(xn - myIntersectionPoints(0, 0),
                                             yn - myIntersectionPoints(1, 0),
                                             zn - myIntersectionPoints(2, 0))(reactant()->species(), neighbor->species()),
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
                                                     Z + Site::nNeighborsLimit(),
                                                     reactant()->species(),
                                                     particle->species());

}

double DiffusionReaction::getSaddleEnergyContributionFromNeighborAt(const uint i, const uint j, const uint k, const uint s1, const uint s2)
{
    return m_saddlePotential(m_saddleFieldIndices[0], m_saddleFieldIndices[1], m_saddleFieldIndices[2])
            (i - m_neighborSetIntersectionPoints(m_saddleFieldIndices[0], m_saddleFieldIndices[1], m_saddleFieldIndices[2])(0, 0),
            j - m_neighborSetIntersectionPoints(m_saddleFieldIndices[0], m_saddleFieldIndices[1], m_saddleFieldIndices[2])(1, 0),
            k - m_neighborSetIntersectionPoints(m_saddleFieldIndices[0], m_saddleFieldIndices[1], m_saddleFieldIndices[2])(2, 0))(s1, s2);
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
        s << Site::info(xD(), yD(), zD(), -m_path[0], -m_path[1], -m_path[2], "O");
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


mat        DiffusionReaction::m_rPowers   = {1.0};
mat        DiffusionReaction::m_strengths = {1.0};

double     DiffusionReaction::m_betaChangeScaleFactor;

field<mat> DiffusionReaction::m_potential;

field<field<mat>>
DiffusionReaction::m_saddlePotential;

field<umat::fixed<3, 2> >
DiffusionReaction::m_neighborSetIntersectionPoints;
