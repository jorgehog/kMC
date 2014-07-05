#pragma once

#include <sys/types.h>

#include <functional>

#include <kMC> //TMP

using namespace std;

namespace WLMC
{

class WLMCWindow;

class WLMCSystem
{
public:
    WLMCSystem(const uint nParticles,
               const uint NX,
               const uint NY,
               const uint NZ,
               const uint movesPerSampling,
               const double flatnessCriterion,
               const uint overlap,
               const uint minWindowSize,
               const uint windowIncrementSize,
               const double *f,
               function<double()> URNG);

    virtual bool isOccupiedLoction(const uint x, const uint y, const uint z) const = 0;

    virtual double getValue(const uint particleIndex) const = 0;

    virtual double getValueDifference(const uint particleIndex, const uint xd, const uint yd, const uint zd) const = 0;

    virtual double getTotalValue() const = 0;

    virtual void changePosition(const uint particleIndex, const uint xd, const uint yd, const uint zd) = 0;


    void sampleWindow(WLMCWindow *window);

    bool doSingleMove(WLMCWindow *window);

    void findDestination(const uint destination, uint &xd, uint &yd, uint &zd);

    void locateGlobalExtremaValues(double &min, double &max, kMC::KMCSolver *solver);

    const uint &nParticles() const
    {
        return m_nParticles;
    }

    const uint &movesPerWindowCheck() const
    {
        return m_movesPerSampling;
    }

    const double &flatnessCriterion() const
    {
        return m_flatnessCriterion;
    }

    const uint &overlap() const
    {
        return m_overlap;
    }

    const uint &minWindowSize() const
    {
        return m_minWindowSize;
    }

    const uint &windowIncrementSize() const
    {
        return m_windowIncrementSize;
    }

    const double &f() const
    {
        return *m_f;
    }

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

    const uint m_movesPerSampling;

    const double m_flatnessCriterion;
    const uint m_overlap;
    const uint m_minWindowSize;
    const uint m_windowIncrementSize;

    const double *m_f;

    const function<double()> m_URNG;


    double getGlobalExtremum(const extrema type);

    void getRandomParticleAndDestination(uint &particleIndex, uint &xd, uint &yd, uint &zd);

    void randomizeParticlePositions();

};

}
