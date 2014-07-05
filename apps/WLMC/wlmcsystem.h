#pragma once

#include <sys/types.h>

#include <functional>

#include <kMC>

using namespace std;

namespace WLMC
{

class WLMCWindow;

class WLMCSystem
{
public:
    WLMCSystem(const uint nParticles, const uint NX, const uint NY, const uint NZ, function<double()> URNG);

    virtual bool isOccupiedLoction(const uint x, const uint y, const uint z) const = 0;

    virtual double getValue(const uint particleIndex) const = 0;

    virtual double getValueDifference(const uint particleIndex, const uint xd, const uint yd, const uint zd) const = 0;

    virtual double getTotalValue() const = 0;


    virtual void clearParticles() = 0;

    virtual void addParticle(const uint x, const uint y, const uint z) = 0;

    virtual void changePosition(const uint particleIndex, const uint xd, const uint yd, const uint zd) = 0;



    void performMove(WLMCWindow *window);

    void findDestination(const uint destination, uint &xd, uint &yd, uint &zd);

    void locateGlobalExtremaValues(double &min, double &max, kMC::KMCSolver *solver);

    enum class extrema
    {
        minimum,
        maximum
    };

private:

    const uint m_nParticles;

    const uint m_NX;
    const uint m_NY;
    const uint m_NZ;

    const uint m_volume;
    const uint m_freeVolume;

    const function<double()> m_URNG;

    double getGlobalExtremum(const extrema type);

    void getRandomParticleAndDestination(uint &particleIndex, uint &xd, uint &yd, uint &zd);

    void randomizeParticlePositions();

};

}
