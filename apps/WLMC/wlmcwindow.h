#pragma once

#include <armadillo>
#include <vector>

#include <kMC> //TMP

using std::vector;

using namespace arma;

namespace WLMC
{

class WLMCSystem;

class WLMCWindow
{
public:

    enum class OVERLAPTYPES
    {
        LOWER,
        UPPER,
        NONE
    };

    WLMCWindow(WLMCSystem *system,
               const vec &parentDOS,
               const vec &parentEnergies,
               const uint lowerLimit,
               const uint upperLimit,
               WLMCWindow::OVERLAPTYPES overlapType);

    WLMCWindow(WLMCSystem *system,
               const uint nBins,
               const double minValue,
               const double maxValue);

    virtual ~WLMCWindow();

    void calculateWindow(kMC::KMCSolver *solver);

    double estimateFlatness(const uint lowerLimit, const uint upperLimit) const;

    void findSubWindows(kMC::KMCSolver *solver);

    uint findFlattestOrigin(const uint lowerLimit, const uint upperLimit) const;

    double getMeanFlatness(const uint lowerLimit, const uint upperLimit) const;

    void findComplementaryRoughAreas(const uint lowerLimitFlat, const uint upperLimitFlat, vector<uvec2> &roughAreas, vector<OVERLAPTYPES> &overlaps) const;

    void findFlatArea(uint &lowerLimit, uint &upperLimit) const;

    void scanForFlattestArea(uint &lowerLimit, uint &upperLimit) const;

    void expandFlattestArea(uint &lowerLimit, uint &upperLimit) const;

    void getSubWindowLimits(WLMCWindow::OVERLAPTYPES overlapType, const uint lowerLimitRough, const uint upperLimitRough, uint &lowerLimit, uint &upperLimit) const;

    void registerVisit(const uint bin);

    uint getBin(double value) const;

    void reset();

    bool isLegal(const double value) const
    {
        return value >= m_minValue && value <= m_maxValue;
    }

    const uint &lowerLimit() const
    {
        return m_lowerLimitOnParent;
    }

    const uint &upperLimit() const
    {
        return m_upperLimitOnParent;
    }

    const vec &DOS() const
    {
        return m_DOS;
    }

    const double &DOS(const uint i) const
    {
        return m_DOS(i);
    }

    const vec &energies() const
    {
        return m_energies;
    }

    const uvec &visitCounts() const
    {
        return m_visitCounts;
    }

    bool isUnsetCount(const uint i) const
    {
        return m_visitCounts(i) == m_unsetCount;
    }

    bool isFlat(const uint lowerLimit, const uint upperLimit) const;

    bool isFlat() const
    {
        return isFlat(0, m_nbins);
    }

    const WLMCWindow::OVERLAPTYPES &overlapType()
    {
        return m_overlapType;
    }

    static constexpr uint m_unsetCount = std::numeric_limits<uint>::max();

private:

    WLMCSystem *m_system;

    vector<WLMCWindow*> m_subWindows;

    const WLMCWindow::OVERLAPTYPES m_overlapType;

    const uint m_lowerLimitOnParent;
    const uint m_upperLimitOnParent;
    const uint m_nbins;

    const double m_minValue;
    const double m_maxValue;
    const double m_valueSpan;

    vec m_DOS;
    const vec m_energies;
    uvec m_visitCounts;

    WLMCWindow *m_lowerNeighbor;
    WLMCWindow *m_upperNeighbor;

    void setNeighbors(WLMCWindow *lowerNeighbor, WLMCWindow *upperNeighbor)
    {
        m_lowerNeighbor = lowerNeighbor;
        m_upperNeighbor = upperNeighbor;
    }

    void mergeWith(WLMCWindow *other);

    uint topIncrement(const uint upperLimit) const;

    uint bottomIncrement(const uint lowerLimit) const;

    void tmp_output(kMC::KMCSolver *solver, const uint lowerLimitFlat, const uint upperLimitFlat, const vector<uvec2> &roughAreas) const; //TMP

};

}
