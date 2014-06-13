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

    const double m_fBegin;
    const double m_fEnd;
    const double m_flatnessCriteria;

    const function<double(double)> m_fIteratorFunction = [] (double fPrev) {return pow(fPrev, 0.75);};


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

    double estimateFlatness() const;

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
