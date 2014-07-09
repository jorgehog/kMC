#pragma once

#include <kMC>

namespace kMC
{

class WLMCEvent : public KMCEvent
{
public:

    WLMCEvent(const uint adaptiveWindows,
              const uint nbins,
              const uint movesPerSampling,
              const double flatnessCritera,
              const uint overlap,
              const uint nbinsOverMinWindowSizeFlat,
              const uint minWindowSizeRough,
              const uint windowIncrementSize,
              const double fBegin,
              const double fEnd) :
        KMCEvent("kMC::WLMC", "", true),
        m_adaptiveWindows(adaptiveWindows),
        m_nbins(nbins),
        m_movesPerWindowCheck(movesPerSampling),
        m_nbinsOverMinWindowSizeFlat(nbinsOverMinWindowSizeFlat),
        m_minWindowSizeRough(minWindowSizeRough),
        m_windowOverlap(overlap),
        m_windowIncrementSize(windowIncrementSize),
        m_fBegin(fBegin),
        m_fEnd(fEnd),
        m_flatnessCriteria(flatnessCritera)
    {

    }

    void initialize();

protected:

    void execute();

private:

    const uint m_adaptiveWindows;

    const uint m_nbins;
    const uint m_movesPerWindowCheck;

    const uint m_nbinsOverMinWindowSizeFlat;
    const uint m_minWindowSizeRough;
    const uint m_windowOverlap;
    const uint m_windowIncrementSize;

    const double m_fBegin;
    const double m_fEnd;
    const double m_flatnessCriteria;

    const function<double(double)> m_fIteratorFunction = [] (double fPrev) {return sqrt(fPrev);};

    vec m_maxEnergies;
    vec m_minEnergies;

    uvec m_visitCounts;
    vec m_DOS;

    double m_f;

    double m_energySpan;

    uint m_nCount;

    uint m_maxBin;
    uint m_minBin;

    uint m_windowMaxBin;
    uint m_windowMinBin;


    void calculateWindow(const uint startBin, const uint endBin);

    uint pointClosestToMean() const;

    void output() const;

    double estimateFlatness(const uvec &visitCounts) const;

    bool isFlat(const uvec &visitCounts) const;

    uint topIncrement(const uint upperLimit, const uint top) const
    {
        uint newTop = upperLimit + m_windowIncrementSize;

        if (newTop > top)
        {
            return top;
        }

        else
        {
            return newTop;
        }

    }

    uint bottomIncrement(const uint lowerLimit, const uint bottom) const
    {
        if (lowerLimit - bottom < m_windowIncrementSize)
        {
            return bottom;
        }

        else
        {
            return lowerLimit - m_windowIncrementSize;
        }

    }

    void findFlatWindow(const uvec &visitCounts, uint &lowerLimit, uint& upperLimit, const uint origin) const;

    void resetCounts()
    {
        m_visitCounts.fill(m_unsetCount);
    }

    void initializeNewCycle();

    void initRandomConfiguration();

    void moveParticle();

    void findAvailableSite(const uint where, uint &xd, uint &yd, uint &zd) const;

    void performTrialMove(SoluteParticle *particle, const uint xd, const uint yd, const uint zd);

    void registerVisit(const uint bin);

    uint getBin(const double energy);

    void prepNextOccupancyLevel();

    void parseForExtrema(const int sign);

    double getEnergyDifference(const SoluteParticle *particle, const uint xd, const uint yd, const uint zd) const;

    void findAllExtrema();

    void findNewEnergyExtrema();

    bool isLegalBin(const uint bin);

    const uint m_nSkipped = 2;
    const uint m_nStart = 32;

    static constexpr uint m_unsetCount = std::numeric_limits<uint>::max();
};


}
