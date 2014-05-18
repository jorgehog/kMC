#include "diffusionreaction.h"

#include "../../kmcsolver.h"

#include "../../boundary/boundary.h"

#include "../../soluteparticle.h"

#include "../../debugger/debugger.h"

#include "../../potential/potential.h"


using namespace kMC;

DiffusionReaction::DiffusionReaction(SoluteParticle *reactant, int dx, int dy, int dz) :
    Reaction(reactant),
    m_lastUsedEsp(UNSET_ENERGY)
{

    m_path[0] = dx;
    m_path[1] = dy;
    m_path[2] = dz;

    m_pathLength = std::round(std::sqrt(dx*dx + dy*dy + dz*dz));

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

    imat::fixed<3, 2> overlapBox;

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


    for (int xPath = -1; xPath <= 1; ++xPath)
    {
        for (int yPath = -1; yPath <= 1; ++yPath)
        {
            for (int zPath = -1; zPath <= 1; ++zPath)
            {
                if ((xPath == 0) && (yPath == 0) && (zPath == 0))
                {
                    continue;
                }

                i = xPath + 1;
                j = yPath + 1;
                k = zPath + 1;

                overlapBox = m_neighborSetIntersectionPoints(i, j, k);

                m_saddlePotential(i, j, k).set_size(overlapBox(0, 1) - overlapBox(0, 0) + 1,
                                                    overlapBox(1, 1) - overlapBox(1, 0) + 1,
                                                    overlapBox(2, 1) - overlapBox(2, 0) + 1);

                for (int dxn = overlapBox(0, 0); dxn <= overlapBox(0, 1); ++dxn)
                {

                    dx = xPath/2.0 - dxn;

                    for (int dyn = overlapBox(1, 0); dyn <= overlapBox(1, 1); ++dyn)
                    {

                        dy = yPath/2.0 - dyn;

                        for (int dzn = overlapBox(2, 0); dzn <= overlapBox(2, 1); ++dzn)
                        {

                            if (dxn == dyn && dyn == dzn && dzn == 0)
                            {
                                continue;
                            }


                            m_saddlePotential(i, j, k)(dxn - overlapBox(0, 0),
                                                       dyn - overlapBox(1, 0),
                                                       dzn - overlapBox(2, 0)).set_size(SoluteParticle::nSpecies(), SoluteParticle::nSpecies());

                            dz = zPath/2.0 - dzn;

                            r2 = dx*dx + dy*dy + dz*dz;

                            for (uint speciesA = 0; speciesA < SoluteParticle::nSpecies(); ++speciesA)
                            {
                                for (uint speciesB = 0; speciesB < SoluteParticle::nSpecies(); ++speciesB)
                                {

                                    rPower = m_rPowers(speciesA, speciesB);
                                    strength = m_strengths(speciesA, speciesB);

                                    m_saddlePotential(i, j, k)(dxn - overlapBox(0, 0),
                                                               dyn - overlapBox(1, 0),
                                                               dzn - overlapBox(2, 0))(speciesA, speciesB) = strength/pow(r2, rPower/2);
                                }
                            }

                        }
                    }
                }


            }
        }
    }

}

double DiffusionReaction::saddlePotential(SoluteParticle *particle)
{
    int X, Y, Z;

    reactant()->distanceTo(particle, X, Y, Z);

    return saddlePotential(m_path[0] + 1, m_path[1] + 1, m_path[2] + 1, X, Y, Z, reactant()->species(), particle->species());
}

double DiffusionReaction::saddlePotential(const uint i,
                                          const uint j,
                                          const uint k,
                                          const int dx,
                                          const int dy,
                                          const int dz,
                                          const uint speciesA,
                                          const uint speciesB)
{
    return getSaddlePot(i, j, k)
            (dx - m_neighborSetIntersectionPoints(i, j, k)(0, 0),
             dy - m_neighborSetIntersectionPoints(i, j, k)(1, 0),
             dz - m_neighborSetIntersectionPoints(i, j, k)(2, 0))
            (speciesA, speciesB);
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

    double Esp = 0;

    //tmp
    if (SoluteParticle::ss != NULL)
    {
        Esp += SoluteParticle::ss->evaluateSaddleFor(this);
    }

    Site *targetSite;
    SoluteParticle *neighbor;


    const field<mat> & saddlePot = getSaddlePot(m_saddleFieldIndices[0],
            m_saddleFieldIndices[1],
            m_saddleFieldIndices[2]);

    const imat::fixed<3, 2> & myIntersectionPoints = m_neighborSetIntersectionPoints(m_saddleFieldIndices[0],
            m_saddleFieldIndices[1],
            m_saddleFieldIndices[2]);

    for (int dxn = myIntersectionPoints(0, 0); dxn <= myIntersectionPoints(0, 1); ++dxn)
    {
        for (int dyn = myIntersectionPoints(1, 0); dyn <= myIntersectionPoints(1, 1); ++dyn)
        {
            for (int dzn = myIntersectionPoints(2, 0); dzn <= myIntersectionPoints(2, 1); ++dzn)
            {

                if (dxn == dyn && dyn == dzn && dzn == 0)
                {
                    continue;
                }

                targetSite = Site::neighborhood(x(), y(), z(), dxn, dyn, dzn);

                if (targetSite == NULL)
                {
                    continue;
                }

                else if (!targetSite->isActive())
                {
                    continue;
                }

                neighbor = targetSite->associatedParticle();

                double dEsp = saddlePot(dxn - myIntersectionPoints(0, 0),
                                        dyn - myIntersectionPoints(1, 0),
                                        dzn - myIntersectionPoints(2, 0))(reactant()->species(), neighbor->species());


                Esp += dEsp;

                KMCDebugger_AssertClose(dEsp, getSaddleEnergyContributionFrom(targetSite->associatedParticle()), 1E-10,
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

    return getSaddleEnergyContributionFromNeighborAt(X, Y, Z, reactant()->species(), particle->species());

}

double DiffusionReaction::getSaddleEnergyContributionFromNeighborAt(const int dxn, const int dyn, const int dzn, const uint s1, const uint s2)
{

    double dx = path(0)/2.0;
    double dy = path(1)/2.0;
    double dz = path(2)/2.0;

    double dxf = dx - dxn;
    double dyf = dy - dyn;
    double dzf = dz - dzn;

    double dr2 = dxf*dxf + dyf*dyf + dzf*dzf;

    double dEsp2 = m_strengths(s1, s2)/std::pow(dr2, m_rPowers(s1, s2)/2);

    return dEsp2;
}


imat::fixed<3, 2> DiffusionReaction::makeSaddleOverlapMatrix(const ivec & relCoor)
{

    imat::fixed<3, 2> overlap;
    int span = static_cast<int>(Site::nNeighborsLimit());

    for (uint xyz = 0; xyz < 3; ++xyz)
    {
        if (relCoor(xyz) == 1)
        {

            overlap(xyz, 0) = -span + 1;
            overlap(xyz, 1) = span;

        }

        else if (relCoor(xyz) == -1)
        {

            overlap(xyz, 0) = -span;
            overlap(xyz, 1) = span - 1;
        }

        else
        {

            overlap(xyz, 0) = -span;
            overlap(xyz, 1) = span;

            KMCDebugger_Assert(relCoor(xyz), ==, 0, "There should be no other option.");
        }
    }

    return overlap;
}

void DiffusionReaction::calcRate()
{

    double newRate;

    KMCDebugger_Assert(updateFlag(), !=, UNSET_UPDATE_FLAG);

    if (updateFlag() == defaultUpdateFlag)
    {

        double Esp = getSaddleEnergy();
        const double &E   = reactant()->energy();

        double pathLength = std::sqrt(m_path[0]*m_path[0] + m_path[1]*m_path[1] + m_path[2]*m_path[2]);

        newRate = linearRateScale()*std::exp(beta()*(E - Esp))/pathLength;

        m_lastUsedEsp = Esp;
    }

    else if(updateFlag() == updateKeepSaddle)
    {
        newRate = rate()*exp(reactant()->energy() - lastUsedEnergy());
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

field<imat::fixed<3, 2> >
DiffusionReaction::m_neighborSetIntersectionPoints;
