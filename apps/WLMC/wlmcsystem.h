#pragma once

#include <sys/types.h>

#include <functional>

using namespace std;

namespace WLMC
{

class WLMCWindow;

class WLMCSystem
{
public:
    WLMCSystem(const uint nParticles, const uint NX, const uint NY, const uint NZ, function<double()> URNG);

    virtual void getLimits(uint &lowerLimit, uint &upperLimit) = 0;

    virtual bool isOccupiedLoction(const uint x, const uint y, const uint z) const = 0;

    virtual double getTotalEnergy() const = 0;

    virtual double getEnergyDifference(const uint particle, const uint xd, const uint yd, const uint zd) const = 0;

    virtual void changePosition(const uint xd, const uint yd, const uint zd) = 0;


    void performMove(WLMCWindow *window);

    void findDestination(const uint destination, uint &xd, uint &yd, uint &zd);

private:

    const uint m_nParticles;

    const uint m_NX;
    const uint m_NY;
    const uint m_NZ;

    const uint m_volume;
    const uint m_freeVolume;

    const function<double()> m_URNG;

};

}
