#pragma once

#include <kMC>

namespace kMC
{

class WLMCEvent : public KMCEvent
{
public:

    WLMCEvent(const uint nbins,
              const double fBegin = datum::e,
              const double fEnd=1.000001,
              const double flatnessCritera = 0.85) :
        KMCEvent("kMC::WLMC", "", true),
        m_nbins(nbins),
        m_minWindow(nbins/10),
        m_windowIncrementSize(5),
        m_fBegin(fBegin),
        m_fEnd(fEnd),
        m_flatnessCriteria(flatnessCritera)
    {

    }

    void initialize();

protected:

    void execute();

private:

    const uint m_nbins;

    const uint m_minWindow;

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
