#pragma once

#include <armadillo>
#include <vector>

using std::vector;

using namespace arma;

namespace WLMC
{

class WLMCSystem;

class WLMCWindow
{
public:
    WLMCWindow(const vec &DOS,
               const uint lowerLimit,
               const uint upperLimit,
               const double minValue,
               const double maxValue);

    void calculateWindow(double f);

    double estimateFlatness() const;

    void findFlatAreas(vector<uvec2> &flatAreas) const;

    void registerVisit(const uint bin);

    uint getBin(double value) const;

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

    static constexpr uint m_unsetCount = std::numeric_limits<uint>::max();

private:

    WLMCSystem *system;

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

    double m_f;

    void performMove();

    void setNeighbors(WLMCWindow *lowerNeighbor, WLMCWindow *upperNeighbor)
    {
        m_lowerNeighbor = lowerNeighbor;
        m_upperNeighbor = upperNeighbor;
    }

    void mergeWith(WLMCWindow *other);

};

}
