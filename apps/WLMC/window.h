#pragma once

#include <armadillo>
#include <vector>

using std::vector;

using namespace arma;

namespace WLMC
{

class System;

struct WindowParams;

class Window
{
public:

    enum class OverlapTypes
    {
        Lower,
        Upper,
        None
    };

    Window(System *system,
           const vec &parentDOS,
           const vec &parentEnergies,
           const uint lowerLimit,
           const uint upperLimit,
           const uint overlapPoint,
           Window::OverlapTypes overlapType, bool allowSubwindowing);

    Window(System *system,
           const uint nBins,
           const double minValue,
           const double maxValue, bool allowSubwindowing);

    Window(const Window &parentWindow, const WindowParams &windowParams);

    virtual ~Window();

    void adapt(const uint lowerLimit, const uint upperLimit);

    vec getHistogram(const uint nmoves);

    void loadInitialConfig();

    void calculateWindow();

    double estimateFlatness(const uint lowerLimit, const uint upperLimit) const;

    void findSubWindows();

    double getMeanFlatness(const uint lowerLimit, const uint upperLimit) const;

    double getMeanFlatness() const
    {
        return getMeanFlatness(0, m_nbins);
    }

    void findComplementaryRoughAreas(vector<WindowParams> &roughWindowParams) const;

    bool findFlatArea();

    bool scanForFlattestArea();

    void expandFlattestArea();

    bool flatProfileIsContinousOnParent() const;


    void getSubWindowLimits(WindowParams &windowParams) const;

    void registerVisit(const uint bin);

    uint getBin(double value) const;

    void reset();

    bool allowsSubwindowing() const;

    bool flatAreaIsInsufficient() const;

    bool flatAreaIsDominant() const;

    bool isLegal(const double value) const
    {
        return value >= m_minValue && value <= m_maxValue;
    }

    const double &minValue() const
    {
        return m_minValue;
    }

    const double &maxValue() const
    {
        return m_maxValue;
    }

    const double &valueSpan() const
    {
        return m_valueSpan;
    }

    const uint &lowerLimitOnParent() const
    {
        return m_lowerLimitOnParent;
    }

    const uint &upperLimitOnParent() const
    {
        return m_upperLimitOnParent;
    }

    System *system() const
    {
        return m_system;
    }

    const vec &DOS() const
    {
        return m_DOS;
    }

    const double &DOS(const uint i) const
    {
        return m_DOS(i);
    }

    void DOS(const vec newDOS)
    {
        m_DOS = newDOS;
    }

    const vec &energies() const
    {
        return m_energies;
    }

    const uvec &visitCounts() const
    {
        return m_visitCounts;
    }

    const uint &visitCounts(const uint i) const
    {
        return m_visitCounts(i);
    }

    const uint &nbins() const
    {
        return m_nbins;
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

    const Window::OverlapTypes &overlapType()
    {
        return m_overlapType;
    }

    static constexpr uint m_unsetCount = std::numeric_limits<uint>::max();

private:

    System *m_system;

    vector<Window*> m_subWindows;
    bool m_allowsSubwindowing;

    const uint m_overlapPoint;
    const Window::OverlapTypes m_overlapType;

    uint m_lowerLimitOnParent;
    uint m_upperLimitOnParent;
    uint m_nbins;

    double m_minValue;
    double m_maxValue;
    double m_valueSpan;

    vec m_DOS;
    vec m_energies;
    uvec m_visitCounts;

    uint m_flatAreaLower;
    uint m_flatAreaUpper;

    void mergeWith(Window *other);

    uint topIncrement(const uint upperLimit) const;

    uint bottomIncrement(const uint lowerLimit) const;

    void tmp_output(const uint lowerLimitFlat, const uint upperLimitFlat, const vector<WindowParams> &roughAreas) const;

};

}
