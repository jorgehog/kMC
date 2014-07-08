#pragma once

#include "wlmcsystem.h"

using namespace WLMC;

namespace kMC
{

class KMCSolver;

class KMCWLMCSystem : public WLMCSystem
{
public:
    KMCWLMCSystem(kMC::KMCSolver *solver,
                  const uint movesPerSampling,
                  const double flatnessCriterion,
                  const uint overlap,
                  const uint nbinsOverMinWindowSizeFlat,
                  const uint minWindowSizeRough,
                  const uint windowIncrementSize,
                  const double *f);

    bool isOccupiedLoction(const uint x, const uint y, const uint z) const;

    double getValue(const uint particleIndex) const;

    double getValueDifference(const uint particleIndex, const uint xd, const uint yd, const uint zd) const;

    double getTotalValue() const;

    void changePosition(const uint particleIndex, const uint xd, const uint yd, const uint zd);

    void getPosition(const uint particleIndex, uint &x, uint &y, uint &z) const;

private:

    KMCSolver *m_solver;

};

}
