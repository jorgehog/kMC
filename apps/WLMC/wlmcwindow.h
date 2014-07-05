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
    WLMCWindow(WLMCSystem *system,
               const vec &DOS,
               const uint lowerLimit,
               const uint upperLimit,
               const double minValue,
               const double maxValue);

    WLMCWindow(WLMCSystem *system,
               const uint nBins,
               const double minValue,
               const double maxValue);

    virtual ~WLMCWindow();

    void calculateWindow(kMC::KMCSolver *solver);

    double estimateFlatness(const uvec &visitCounts) const;

    void findFlatAreas(vector<uvec2> &flatAreas, const uint lowerLimit, const uint upperLimit) const;

    uint findFlattestOrigin(const uint lowerLimit, const uint upperLimit) const;

    double getMeanFlatness(const uint lowerLimit, const uint upperLimit) const;

    void findComplementaryRoughAreas(const vector<uvec2> &flatAreas, vector<uvec2> &roughAreas, const uint lowerLimit, const uint _upperLimit) const;

    void findFlatArea(uint &upperLimit, uint &lowerLimit, const uint origin) const;

    void registerVisit(const uint bin);

    uint getBin(double value) const;

    void reset();

    bool isLegal(const uint bin) const
    {
        return bin >= m_lowerLimit && bin <= m_upperLimit;
    }

    const double &DOS(const uint i) const
    {
        return m_DOS(i);
    }

    bool isUnsetCount(const uint i) const
    {
        return m_visitCounts(i) == m_unsetCount;
    }

    bool isFlat(const uvec& visitCounts) const;

    bool isFlat() const
    {
        return isFlat(m_visitCounts);
    }

    static constexpr uint m_unsetCount = std::numeric_limits<uint>::max();

private:

    WLMCSystem *m_system;

    const uint m_lowerLimit;
    const uint m_upperLimit;
    const uint m_nbins;

    const double m_minValue;
    const double m_maxValue;
    const double m_valueSpan;

    vec m_DOS;
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

    void tmp_output(kMC::KMCSolver *solver, const vector<uvec2> &flatAreas, const vector<uvec2> &roughAreas) const; //TMP

};

}
