#pragma once

#include <armadillo>

#include <functional>

using namespace arma;
using namespace std;

namespace WLMC
{

class Window;

class System
{
public:
    System(const uint nParticles,
               const uint NX,
               const uint NY,
               const uint NZ,
               const uint movesPerSampling,
               const double flatnessCriterion,
               const uint overlap,
               const uint nbinsOverMinWindowSizeFlat,
               const uint minWindowSize,
               const uint windowIncrementSize,
               const double *f,
               function<double()> URNG);

    virtual bool isOccupiedLoction(const uint x, const uint y, const uint z) const = 0;

    virtual double getValue(const uint particleIndex) const = 0;

    virtual double getValueDifference(const uint particleIndex, const uint xd, const uint yd, const uint zd) const = 0;

    virtual double getTotalValue() const = 0;

    virtual void changePosition(const uint particleIndex, const uint xd, const uint yd, const uint zd) = 0;

    virtual void getPosition(const uint particleIndex, uint &x, uint &y, uint &z) const = 0;

    void sampleWindow(Window *window);

    bool doWLMCMove(Window *window);

    void doRandomMove();

    void findDestination(const uint destination, uint &xd, uint &yd, uint &zd) const;

    void locateGlobalExtremaValues(double &min, double &max);

    void setupPresetWindowConfigurations(const double min, const double max, const uint n);

    void loadConfigurationClosestToValue(const double value);

    uint getPresetBinFromValue(const double value) const;

    void clipWindow(Window &window) const;

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

    uint minWindowSizeFlat(const uint nbins) const
    {
        return nbins/m_nbinsOverMinWindowSizeFlat;
    }

    uint minWindowSize() const
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
    const uint m_nbinsOverMinWindowSizeFlat;
    const uint m_minWindowSize;
    const uint m_windowIncrementSize;

    const double *m_f;

    const function<double()> m_URNG;

    ucube m_presetWindowConfigurations;
    vec m_presetWindowValues;

    double getGlobalExtremum(const extrema type);

    void getRandomParticleAndDestination(uint &particleIndex, uint &xd, uint &yd, uint &zd) const;

    void randomizeParticlePositions();

};

}
