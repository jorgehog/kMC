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
    WLMCWindow(WLMCSystem *system,
               const vec &DOS,
               const uint lowerLimit,
               const uint upperLimit,
               const double *f,
               const double minValue,
               const double maxValue);

    WLMCWindow(WLMCSystem *system,
               const uint nBins,
               const double *f,
               const double minValue,
               const double maxValue);

    virtual ~WLMCWindow();

    void calculateWindow();

    double estimateFlatness() const;

    void findFlatAreas(vector<uvec2> &flatAreas) const;

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

    static constexpr uint m_unsetCount = std::numeric_limits<uint>::max();

private:

    WLMCSystem *m_system;

    const double* m_f;

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


    const double &f() const
    {
        return *m_f;
    }

    void setNeighbors(WLMCWindow *lowerNeighbor, WLMCWindow *upperNeighbor)
    {
        m_lowerNeighbor = lowerNeighbor;
        m_upperNeighbor = upperNeighbor;
    }

    void mergeWith(WLMCWindow *other);

};

}
