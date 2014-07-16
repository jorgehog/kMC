#include "window.h"
#include "system.h"
#include "windowparams.h"

#include <kMC>
#include <BADAss/badass.h>

using namespace WLMC;

Window::Window(System *system,
               const vec &parentDOS,
               const vec &parentEnergies,
               const uint lowerLimit,
               const uint upperLimit,
               const uint overlapPoint,
               OverlapTypes overlapType,
               bool allowSubwindowing) :
    m_system(system),
    m_allowsSubwindowing(allowSubwindowing),
    m_overlapPoint(overlapPoint),
    m_overlapType(overlapType),
    m_lowerLimitOnParent(lowerLimit),
    m_upperLimitOnParent(upperLimit),
    m_nbins(upperLimit - lowerLimit),
    m_minValue(parentEnergies(lowerLimit)),
    m_maxValue(parentEnergies(upperLimit - 1)),
    m_valueSpan(m_maxValue - m_minValue),
    m_DOS(parentDOS(span(lowerLimit, upperLimit - 1))),
    m_energies(parentEnergies(span(lowerLimit, upperLimit - 1))),
    m_visitCounts(uvec(m_nbins))
{
    m_visitCounts.fill(m_unsetCount);

    BADAss(upperLimit, >, lowerLimit, "illegal window.");
}

Window::Window(System *system, const uint nBins,
               const double minValue,
               const double maxValue,
               bool allowSubwindowing) :
    Window(system,
           ones(nBins),
           linspace(minValue, maxValue, nBins),
           0,
           nBins,
           0,
           Window::OverlapTypes::None,
           allowSubwindowing)
{

}

Window::Window(const Window &parentWindow, const WindowParams &windowParams) :
    Window(parentWindow.system(),
           parentWindow.DOS(),
           parentWindow.energies(),
           windowParams.m_lowerLimit,
           windowParams.m_upperLimit,
           windowParams.m_overlapPoint,
           windowParams.m_overlapType,
           windowParams.m_allowSubwindowing)
{

}

Window::~Window()
{
    for (Window *window : m_subWindows)
    {
        delete window;
    }

    m_subWindows.clear();

    m_DOS.clear();
    m_visitCounts.clear();
    m_system = NULL;
}

void Window::adapt(const uint lowerLimit, const uint upperLimit)
{
    m_DOS = m_DOS(span(lowerLimit, upperLimit - 1));
    m_energies = m_energies(span(lowerLimit, upperLimit - 1));
    m_visitCounts = m_visitCounts(span(lowerLimit, upperLimit - 1));

    m_nbins = upperLimit - lowerLimit;

    m_minValue = m_energies(0);
    m_maxValue = m_energies(m_nbins - 1);

    m_valueSpan = m_maxValue - m_minValue;

    m_lowerLimitOnParent += lowerLimit;
    m_upperLimitOnParent -= m_upperLimitOnParent - upperLimit;
}

vec Window::getHistogram(const uint nmoves)
{
    double value = 0;

    vec histogram(m_nbins, fill::zeros);

    uint move = 0;
    while (move < nmoves)
    {
        m_system->doRandomMove();

        value = m_system->getTotalValue();

        if (isLegal(value))
        {
            uint bin = getBin(m_system->getTotalValue());

            histogram(bin)++;
        }

        move++;
    }

    return normalise(histogram);

}

void Window::loadInitialConfig()
{
    m_system->loadConfigurationClosestToValue((m_maxValue + m_minValue)/2);
}

void Window::calculateWindow()
{

    cout << "sampling on " << m_lowerLimitOnParent << " " << m_upperLimitOnParent << " f = " << m_system->f() << endl;

    vector<WindowParams> roughWindowParams;

    while (!isFlat())
    {

        m_system->sampleWindow(this);

        m_DOS = normalise(m_DOS);

        findSubWindows();

        for (Window *subWindow : m_subWindows)
        {
            subWindow->loadInitialConfig();
            subWindow->calculateWindow();
            mergeWith(subWindow);
        }


        BADAssBool(!m_subWindows.empty() && isFlat(), "subwindows were not merged properly.");

    }

}

double Window::estimateFlatness(const uint lowerLimit, const uint upperLimit) const
{

    BADAss(upperLimit - lowerLimit, <, m_system->minWindowSize(),
           "Flatness requested on window lower than minimum limit of rough areas.");

    uint min = std::numeric_limits<uint>::max();

    double mean = 0;

    uint nSetCounts = 0;

    for (uint i = lowerLimit; i < upperLimit; ++i)
    {
        uint vc = m_visitCounts(i);

        if (!(vc == m_unsetCount))
        {

            if (vc < min)
            {
                min = vc;
            }

            mean += vc;

            nSetCounts++;
        }
    }

    if (nSetCounts == 0)
    {
        return 0;
    }

    mean /= nSetCounts;

    return min/mean;

}

void Window::findSubWindows()
{
    if (!allowsSubwindowing())
    {
        return;
    }

    else if (!findFlatArea())
    {
        return;
    }

    else if (!flatProfileIsContinousOnParent())
    {
        return;
    }

    else if (flatAreaIsDominant())
    {
        m_allowsSubwindowing = false;
        return;
    }

    else if (flatAreaIsInsufficient())
    {
        return;
    }

    else
    {

        vector<WindowParams> roughWindowParams;
        findComplementaryRoughAreas(roughWindowParams);

        for (WindowParams &windowParam : roughWindowParams)
        {
            getSubWindowLimits(windowParam);

            m_subWindows.push_back(new Window(*this, windowParam));

        }

    }




















    //    bool allowSubWindowing;
    //    uint lowerLimitNewWindow, upperLimitNewWindow;

    //    if (isFlat(lowerLimitFlat, upperLimitFlat))
    //    {
    //        vector<uvec2> roughAreas;
    //        vector<Window::OVERLAPTYPES> overlaps;

    //        findComplementaryRoughAreas(lowerLimitFlat, upperLimitFlat, roughAreas, overlaps);

    //        vector<uvec2> _roughAreas(roughAreas.size());
    //        vector<uint> overlapPointBins(roughAreas.size());

    //        for (uint i = 0; i < roughAreas.size(); ++i)
    //        {
    //            uvec2 roughArea = roughAreas.at(i);
    //            Window::OVERLAPTYPES overlapType = overlaps.at(i);

    //            if (roughArea(1) - roughArea(0) < m_system->minWindowSizeRough())
    //            {
    //                cout << "too low rought area" << endl;
    //                return;
    //            }

    //            if (overlapType == Window::OVERLAPTYPES::LOWER)
    //            {
    //                overlapPointBins.at(i) = roughArea(1);
    //            }
    //            else
    //            {
    //                overlapPointBins.at(i) = roughArea(0);
    //            }

    //            getSubWindowLimits(overlapType, roughArea(0), roughArea(1), lowerLimitNewWindow, upperLimitNewWindow, allowSubWindowing);

    //            _roughAreas.at(i) = {lowerLimitNewWindow, upperLimitNewWindow};
    //        }

    //        tmp_output(lowerLimitFlat, upperLimitFlat, _roughAreas);

    //        for (uint i = 0; i < roughAreas.size(); ++i)
    //        {
    //            uvec2 roughArea = roughAreas.at(i);
    //            Window::OVERLAPTYPES overlapType = overlaps.at(i);


    //            getSubWindowLimits(overlapType, roughArea(0), roughArea(1), lowerLimitNewWindow, upperLimitNewWindow, allowSubWindowing);

    //            cout << "############### SPAWNED SUB WINDOW ###################### " <<  lowerLimitNewWindow << " " << upperLimitNewWindow << endl;

    //            m_subWindows.push_back(new Window(m_system, m_DOS, m_energies, lowerLimitNewWindow, upperLimitNewWindow, overlapPointBins.at(i), overlapType, allowSubWindowing));
    //        }

    //        return;
    ////        return getMeanFlatness(lowerLimitFlat, upperLimitFlat);

    //    }
    //    else
    //    {
    //        tmp_output(0, 0, {{0, m_nbins}});
    //        return;
    //    }

}


double Window::getMeanFlatness(const uint lowerLimit, const uint upperLimit) const
{
    double mean = 0;

    uint nSetCounts = 0;

    for (uint i = lowerLimit; i < upperLimit; ++i)
    {
        uint vc = m_visitCounts(i);

        if (!(vc == m_unsetCount))
        {
            mean += vc;

            nSetCounts++;
        }
    }

    return mean/nSetCounts;
}

void Window::findComplementaryRoughAreas(vector<WindowParams> &roughWindowParams) const
{

    if (m_flatAreaLower != 0)
    {
        roughWindowParams.push_back(WindowParams(0, m_flatAreaLower, Window::OverlapTypes::Upper));
    }

    if (m_flatAreaUpper != m_nbins)
    {
        roughWindowParams.push_back(WindowParams(m_flatAreaUpper, m_nbins, Window::OverlapTypes::Lower));
    }

}

bool Window::findFlatArea()
{
    m_flatAreaLower = m_flatAreaUpper = 0;

    if (scanForFlattestArea())
    {
        expandFlattestArea();
        return true;
    }

    return false;
}

bool Window::scanForFlattestArea()
{
    //performs maximization of flatness for windows of size minWindowSize for increments
    //of size windowIncrementSize. Interval of maximum is stored in lowerLimit and upperLimit.

    uint lowerLimitScan = 0;
    uint upperLimitScan = m_system->minWindowSizeFlat(m_nbins);

    double Fn;

    constexpr double unsetFn = -1;

    double maxFn = unsetFn;

    while (upperLimitScan < m_nbins)
    {
        Fn = estimateFlatness(lowerLimitScan, upperLimitScan);

        if (Fn > maxFn)
        {
            maxFn = Fn;

            m_flatAreaLower = lowerLimitScan;
            m_flatAreaUpper = upperLimitScan;

        }

        lowerLimitScan += m_system->windowIncrementSize();
        upperLimitScan += m_system->windowIncrementSize();
    }

    return maxFn < m_system->flatnessCriterion();

}

void Window::expandFlattestArea()
{
    //performs maximization of s = u - l under criteria that [u, l] is flat. originates at [lowerLimit, upperLimit]
    //resulting interval where the max of s is found is stored in lowerLimit and upperLimit.

    int expandingLowerLimit = m_flatAreaLower;
    uint upperLimitStart = m_flatAreaUpper;

    uint maxSpan = m_flatAreaUpper - m_flatAreaLower;

    while (expandingLowerLimit > 0)
    {
        uint expandingUpperLimit = upperLimitStart;

        while (expandingUpperLimit <= m_nbins)
        {
            if (isFlat(expandingLowerLimit, expandingUpperLimit))
            {
                uint span = expandingUpperLimit - expandingLowerLimit;

                if (span > maxSpan)
                {
                    maxSpan = span;

                    m_flatAreaLower = expandingLowerLimit;
                    m_flatAreaUpper = expandingUpperLimit;

                }
            }

            expandingUpperLimit += m_system->windowIncrementSize();
        }

        expandingLowerLimit -= m_system->windowIncrementSize();
    }

}

bool Window::flatProfileIsContinousOnParent() const
{
    if (m_overlapType == Window::OverlapTypes::Lower)
    {
        return m_flatAreaLower == 0;
    }

    else if (m_overlapType == Window::OverlapTypes::Upper)
    {
        return m_flatAreaUpper == m_nbins;
    }

    else
    {
        return true;
    }
}

void Window::getSubWindowLimits(WindowParams &windowParams) const
{









//    uint newSize = upperLimitRough - lowerLimitRough + m_system->overlap();

//    //If the desired area is lower than the minimum window, we need to increase it.
//    if (upperLimitRough - lowerLimitRough < m_system->minWindowSize() || newSize < m_system->minWindowSizeFlat(newSize))
//    {
//        allowSubWindowing = false;
//        return;
//    }

//    else
//    {
//        lowerLimit = lowerLimitRough;
//        upperLimit = upperLimitRough;

//        allowSubWindowing = true;
//    }

//    //We are not guaranteed a window of size at least minWindow.
//    //Now we add the overlap.

//    if (overlapType == Window::OverlapTypes::Lower)
//    {
//        if (lowerLimit < m_system->overlap())
//        {
//            cout << "ERROR IN LOWER OVERLAP " << lowerLimitRough << " " << upperLimitRough << " " << lowerLimit << endl;
//            throw std::runtime_error("asds");
//        }

//        lowerLimit -= m_system->overlap();
//    }
//    else
//    {
//        if (upperLimit > m_nbins - m_system->overlap())
//        {
//            cout << "ERROR IN UPPER OVERLAP " << lowerLimitRough << " " << upperLimitRough << " " << upperLimit << endl;
//            throw std::runtime_error("asds");
//        }

//        upperLimit += m_system->overlap();
//    }

}

void Window::registerVisit(const uint bin)
{
    if (isUnsetCount(bin))
    {
        m_visitCounts(bin) = 1;
    }
    else
    {
        m_visitCounts(bin)++;
    }

    m_DOS(bin) *= m_system->f();
}

uint Window::getBin(double value) const
{
    uint bin = m_nbins*(value - m_minValue)/m_valueSpan;

    if (bin == m_nbins)
    {
        return m_nbins - 1;
    }

    return bin;
}

void Window::reset()
{
    for (Window *subWindow : m_subWindows)
    {
        delete subWindow;
    }

    m_subWindows.clear();

    m_visitCounts.fill(m_unsetCount);
}

bool Window::allowsSubwindowing() const
{
    return m_allowsSubwindowing;
}

bool Window::flatAreaIsInsufficient() const
{
    uint N = 0;
    if (m_overlapType == OverlapTypes::None)
    {
        N = 2;
    }
    else
    {
        N = 1;
    }

    return m_flatAreaUpper - m_flatAreaLower < N*m_system->overlap();
}

bool Window::flatAreaIsDominant() const
{
    return (m_nbins - (m_flatAreaUpper - m_flatAreaLower)) < m_system->minWindowSize();
}

bool Window::isFlat(const uint lowerLimit, const uint upperLimit) const
{
    return estimateFlatness(lowerLimit, upperLimit) >= m_system->flatnessCriterion();
}

void Window::mergeWith(Window *other)
{

    double meanVisitAtFlatArea = getMeanFlatness(m_flatAreaLower, m_flatAreaUpper);

    cout << "############### MERGED SUB WINDOW ###################### " <<  other->lowerLimitOnParent() << " " << other->upperLimitOnParent() << endl;

    uint pointSub;
    uint pointSup;

    uint repl_l;
    uint repl_u;

    if (other->overlapType() == Window::OverlapTypes::Lower)
    {
        pointSub = m_system->overlap();
        repl_l = other->lowerLimitOnParent() + m_system->overlap();
        repl_u = other->upperLimitOnParent();
    }

    else
    {
        pointSub = other->nbins() - m_system->overlap();
        repl_l = other->lowerLimitOnParent();
        repl_u = other->upperLimitOnParent() - m_system->overlap();
    }

    pointSup = other->lowerLimitOnParent() + pointSub;

    double shiftSubSup = m_DOS(pointSup)/other->DOS(pointSub);

    m_DOS(span(repl_l, repl_u - 1)) = shiftSubSup*other->DOS()(span(repl_l - other->lowerLimitOnParent(), repl_u - other->lowerLimitOnParent() - 1));


    cout << "merged " << other->lowerLimitOnParent() << " " << other->upperLimitOnParent() << " with shift on " << shiftSubSup << endl;

    double shiftVisit = meanVisitAtFlatArea/other->getMeanFlatness();

    for (uint i = 0; i < other->nbins(); ++i)
    {
        if (!other->isUnsetCount(i) || isUnsetCount(other->lowerLimitOnParent() + i))
        {
            m_visitCounts(other->lowerLimitOnParent() + i) = other->visitCounts(i)*shiftVisit;
        }
    }

}

uint Window::topIncrement(const uint upperLimit) const
{
    uint newTop = upperLimit + m_system->windowIncrementSize();

    if (newTop > m_nbins)
    {
        return m_nbins;
    }

    else
    {
        return newTop;
    }

}

uint Window::bottomIncrement(const uint lowerLimit) const
{
    if (lowerLimit < m_system->windowIncrementSize())
    {
        return 0;
    }

    else
    {
        return lowerLimit - m_system->windowIncrementSize();
    }

}

void Window::tmp_output(const uint lowerLimitFlat, const uint upperLimitFlat, const vector<WindowParams> &roughAreas) const
{

    ofstream fFlat;
    fFlat.open(kMC::KMCSolver::instance()->filePath() + "/flatness.txt");

    fFlat << lowerLimitFlat << " " << upperLimitFlat << endl;

    fFlat.close();

    ofstream fRough;
    fRough.open(kMC::KMCSolver::instance()->filePath() + "/roughness.txt");
    for (const WindowParams & roughArea : roughAreas)
    {
        fRough << roughArea.m_lowerLimit << " " << roughArea.m_upperLimit << endl;
    }
    fRough.close();


    stringstream file;
    file << kMC::KMCSolver::instance()->filePath() << "/stateDensity" << m_system->nParticles() << ".arma";
    vec E = linspace<vec>(m_minValue, m_maxValue, m_nbins);

    uvec indices = find(m_visitCounts != m_unsetCount);
    uvec vc = m_visitCounts(indices);

    vec dos = m_DOS(indices);
    vec e = E(indices);
    vec idx = conv_to<vec>::from(indices);
    vec vcd = conv_to<vec>::from(vc);

    if (!join_rows(join_rows(join_rows(e, dos), vcd), idx).eval().save(file.str()))
    {
        cout << "failed at storing " << file.str() << endl;
    }
}
