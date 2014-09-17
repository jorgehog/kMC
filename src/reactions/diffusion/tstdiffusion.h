#pragma once

#include "diffusionreaction.h"

namespace kMC
{

class TSTDiffusion : public DiffusionReaction
{
public:

    TSTDiffusion(SoluteParticle *reactant, int dx, int dy, int dz);

    double calcRate();

    void reset();

    void setDirectUpdateFlags(const SoluteParticle *changedReactant, const uint level);

    static void clearAll()
    {
        DiffusionReaction::clearAll();

        m_saddlePotential.reset_objects();
        m_saddlePotential.reset();
        m_neighborSetIntersectionPoints.reset_objects();
        m_neighborSetIntersectionPoints.reset();
    }

    static void setupSaddlePotential();

    double saddlePotential(SoluteParticle *particle);

    double getSaddleEnergy();

    double getSaddleEnergyContributionFrom(const SoluteParticle *particle);

    double getSaddleEnergyContributionFromNeighborAt(const int dxn, const int dyn, const int dzn, const uint s1, const uint s2);

    static imat::fixed<3, 2> makeSaddleOverlapMatrix(const ivec &relCoor);

    static double saddlePotential(const uint i,
                                  const uint j,
                                  const uint k,
                                  const int dx,
                                  const int dy,
                                  const int dz, const uint speciesA, const uint speciesB);

    static const field<mat> & getSaddlePot(const uint i, const uint j, const uint k)
    {
        return m_saddlePotential(i, j, k);
    }

    static const imat::fixed<3, 2> & neighborSetIntersectionPoints(const uint i, const uint j, const uint k)
    {
        return m_neighborSetIntersectionPoints(i, j, k);
    }

    const double & lastUsedEsp() const
    {
        return m_lastUsedEsp;
    }

    enum
    {
        updateKeepSaddle = 1
    };


private:

    static field<field<mat>> m_saddlePotential;

    static field<imat::fixed<3, 2> > m_neighborSetIntersectionPoints;

    double m_lastUsedEsp;

};

}
