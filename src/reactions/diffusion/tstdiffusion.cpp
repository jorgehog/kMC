#include "tstdiffusion.h"

#include "soluteparticle.h"

#include "../../boundary/boundary.h"

#include "../../soluteparticle.h"

#include "../../debugger/debugger.h"

#include "../../potential/potential.h"


using namespace kMC;

TSTDiffusion::TSTDiffusion(SoluteParticle *reactant, int dx, int dy, int dz) :
    DiffusionReaction(reactant, dx, dy, dz),
    m_lastUsedEsp(UNSET_ENERGY)
{

}

void TSTDiffusion::calcRate()
{

    double newRate = 0;
    const double &E   = reactant()->energy();

    BADAss(updateFlag(), !=, UNSET_UPDATE_FLAG);

    if (updateFlag() == defaultUpdateFlag)
    {

        double Esp = getSaddleEnergy();

        newRate = linearRateScale()*std::exp(beta()*(Esp - E))/pathLength();

        m_lastUsedEsp = Esp;
    }

    else if(updateFlag() == updateKeepSaddle)
    {
        BADAss(rate(), !=, UNSET_RATE ,"Saddle can't update when the rate has not been calculated.", KMCBAI( getFinalizingDebugMessage()));
        BADAss(updateFlag(), ==, updateKeepSaddle, "Errorous updateFlag.", KMCBAI( getFinalizingDebugMessage()));
        BADAss(lastUsedEnergy(), !=, UNSET_ENERGY, "energy never calculated before.", KMCBAI( getFinalizingDebugMessage()));

        newRate = rate()*std::exp(beta()*(lastUsedEnergy() - E));

        BADAssClose(getSaddleEnergy(), m_lastUsedEsp, 1E-10, "Saddle energy was not conserved as assumed by flag. ", KMCBAI( getFinalizingDebugMessage()));

    }

    setRate(newRate);

}

void TSTDiffusion::reset()
{

    DiffusionReaction::reset();

    m_lastUsedEsp = UNSET_ENERGY;

}

void TSTDiffusion::setDirectUpdateFlags(const SoluteParticle *changedReactant, const uint level)
{

    BADAss(changedReactant, !=, reactant());

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
                BADAss(Site::nNeighborsLimit() + 1, ==,  d_maxDistance);
                registerUpdateFlag(updateKeepSaddle);
            }

            else
            {
                registerUpdateFlag(defaultUpdateFlag);
            }
        }

    }



}

void TSTDiffusion::setupSaddlePotential()
{
    using std::pow;

    uint i, j, k;
    double rPower, strength;

    double dx, dy, dz, r2;
    imat::fixed<3, 2> overlapBox;

    m_saddlePotential.reset_objects();
    m_saddlePotential.reset();
    m_saddlePotential.set_size(3, 3, 3);

    m_neighborSetIntersectionPoints.reset_objects();
    m_neighborSetIntersectionPoints.reset();
    m_neighborSetIntersectionPoints.set_size(3, 3, 3);


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

                overlapBox = makeSaddleOverlapMatrix({xPath, yPath, zPath});

                m_neighborSetIntersectionPoints(i, j, k) = overlapBox;

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

                                    rPower = DiffusionReaction::rPower(speciesA, speciesB);
                                    strength = DiffusionReaction::strength(speciesA, speciesB);

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

double TSTDiffusion::saddlePotential(SoluteParticle *particle)
{
    int X, Y, Z;

    reactant()->distanceTo(particle, X, Y, Z);

    return saddlePotential(pathIndex(0), pathIndex(1), pathIndex(2), X, Y, Z, reactant()->species(), particle->species());
}

double TSTDiffusion::saddlePotential(const uint i,
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



double TSTDiffusion::getSaddleEnergy()
{

    double Esp = 0;

    //tmp
    if (SoluteParticle::ss != NULL)
    {
        Esp += SoluteParticle::ss->evaluateSaddleFor(this);
    }

    if (reactant()->nNeighborsSum() == 0)
    {
        return Esp;
    }

    Site *targetSite;
    SoluteParticle *neighbor;


    const field<mat> & saddlePot = getSaddlePot(pathIndex(0),
                                                pathIndex(1),
                                                pathIndex(2));

    const imat::fixed<3, 2> & myIntersectionPoints = m_neighborSetIntersectionPoints(pathIndex(0),
                                                                                     pathIndex(1),
                                                                                     pathIndex(2));

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

                BADAssClose(dEsp, getSaddleEnergyContributionFrom(targetSite->associatedParticle()), 1E-10,
                            "Mismatch in saddle energy contribution.",
                            KMCBAI(getFinalizingDebugMessage()));
            }
        }
    }


    return Esp;

}

double TSTDiffusion::getSaddleEnergyContributionFrom(const SoluteParticle *particle)
{
    int X, Y, Z;

    reactant()->distanceTo(particle, X, Y, Z);

    return getSaddleEnergyContributionFromNeighborAt(X, Y, Z, reactant()->species(), particle->species());

}

double TSTDiffusion::getSaddleEnergyContributionFromNeighborAt(const int dxn, const int dyn, const int dzn, const uint s1, const uint s2)
{

    double dx = path(0)/2.0;
    double dy = path(1)/2.0;
    double dz = path(2)/2.0;

    double dxf = dx - dxn;
    double dyf = dy - dyn;
    double dzf = dz - dzn;

    double dr2 = dxf*dxf + dyf*dyf + dzf*dzf;

    double dEsp2 = strength(s1, s2)/std::pow(dr2, rPower(s1, s2)/2);

    return dEsp2;
}


imat::fixed<3, 2> TSTDiffusion::makeSaddleOverlapMatrix(const ivec & relCoor)
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

            BADAss(relCoor(xyz), ==, 0, "There should be no other option.");
        }
    }

    return overlap;
}


field<field<mat> >
TSTDiffusion::m_saddlePotential;

field<imat::fixed<3, 2> >
TSTDiffusion::m_neighborSetIntersectionPoints;
