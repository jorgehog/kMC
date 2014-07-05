#pragma once

#include "wlmcsystem.h"

using namespace WLMC;

namespace kMC
{

class KMCSolver;

class KMCWLMCSystem : public WLMCSystem
{
public:
    KMCWLMCSystem(kMC::KMCSolver *solver);

    bool isOccupiedLoction(const uint x, const uint y, const uint z) const;

    double getValue(const uint particleIndex) const;

    double getValueDifference(const uint particleIndex, const uint xd, const uint yd, const uint zd) const;

    double getTotalValue() const;

    void clearParticles();

    void addParticle(const uint x, const uint y, const uint z);

    void changePosition(const uint particleIndex, const uint xd, const uint yd, const uint zd);

private:

    KMCSolver *m_solver;

};

}
