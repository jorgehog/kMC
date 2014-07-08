#include "wlmcwindow.h"
#include "wlmcsystem.h"

#include <kMC>

using namespace WLMC;

WLMCWindow::WLMCWindow(WLMCSystem *system,
                       const vec &parentDOS,
                       const vec &parentEnergies,
                       const uint lowerLimit,
                       const uint upperLimit,
                       OVERLAPTYPES overlapType) :
    m_system(system),
    m_overlapType(overlapType),
    m_lowerLimitOnParent(lowerLimit),
    m_upperLimitOnParent(upperLimit),
    m_nbins(upperLimit - lowerLimit),
    m_minValue(parentEnergies(lowerLimit)),
    m_maxValue(parentEnergies(m_nbins - 1)),
    m_valueSpan(m_maxValue - m_minValue),
    m_DOS(parentDOS(span(lowerLimit, upperLimit - 1))),
    m_energies(parentEnergies(span(lowerLimit, upperLimit - 1))),
    m_visitCounts(uvec(m_nbins))
{
    m_visitCounts.fill(m_unsetCount);
}

WLMCWindow::WLMCWindow(WLMCSystem *system, const uint nBins,
                       const double minValue,
                       const double maxValue) :
    WLMCWindow(system, ones(nBins), linspace(minValue, maxValue, nBins), 0, nBins, WLMCWindow::OVERLAPTYPES::NONE)
{

}

WLMCWindow::~WLMCWindow()
{
    for (WLMCWindow *window : m_subWindows)
    {
        delete window;
    }

    m_subWindows.clear();

    m_DOS.clear();
    m_visitCounts.clear();
    m_system = NULL;
}

vec WLMCWindow::getHistogram(const uint nmoves)
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

void WLMCWindow::calculateWindow()
{

    cout << "sampling on " << m_lowerLimitOnParent << " " << m_upperLimitOnParent << " f = " << m_system->f() << endl;

    uint lowerLimitFlat = 0;
    uint upperLimitFlat = 0;

    vector<uvec2> roughAreas;
    vector<WLMCWindow::OVERLAPTYPES> _o;


    uint start, end;
    while (!isFlat() || m_subWindows.size() != 0)
    {

        m_system->sampleWindow(this);

        m_DOS = normalise(m_DOS);

        //NO WINDOWS
        //        findFlatArea(lowerLimitFlat, upperLimitFlat);

        //        if (!isFlat(lowerLimitFlat, upperLimitFlat))
        //        {
        //            tmp_output(solver, 0, 0, {{0, m_nbins}});
        //            continue;
        //        }

        //        findComplementaryRoughAreas(lowerLimitFlat, upperLimitFlat, roughAreas, _o);

        //        for (uint i = 0; i < roughAreas.size(); ++i)
        //        {
        //            getSubWindowLimits(_o.at(i), roughAreas.at(i)(0), roughAreas.at(i)(1), start, end);

        //            roughAreas.at(i) = {start, end};
        //        }

        //        tmp_output(solver, lowerLimitFlat, upperLimitFlat, roughAreas);

        //        roughAreas.clear();
        //        _o.clear();

        //WINDOWS


        if (m_system->minWindowSizeRough() > m_system->minWindowSizeFlat(m_nbins))
        {
            cout << "rough window " << m_system->minWindowSizeRough() << " too small for flat window " << m_system->minWindowSizeFlat(m_nbins) << " bins = " << m_nbins << endl;
            tmp_output(0, 0, {{0, m_nbins}});
            continue;
        }

        double mean = findSubWindows();

        for (WLMCWindow *subWindow : m_subWindows)
        {
            subWindow->calculateWindow();
            mergeWith(subWindow, mean);

            findFlatArea(start, end);
            findComplementaryRoughAreas(start, end, roughAreas, _o);
            tmp_output(start, end, roughAreas);

        }

    }

    if (m_subWindows.size() != 0)
    {

    }
    //todo: work out subwindows issues
    //todo: debug subwindows
    //todo: optimize and set transparent parameters for subwindows
    //todo: debug subwindows

}

double WLMCWindow::estimateFlatness(const uint lowerLimit, const uint upperLimit) const
{
    if (upperLimit - lowerLimit < m_system->minWindowSizeRough())
    {
        throw std::runtime_error("error. Flatness requested on window lower than minimum limit of rough areas.");
    }

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

double WLMCWindow::findSubWindows()
{

    uint upperLimitFlat = 0;
    uint lowerLimitFlat = 0;

    findFlatArea(lowerLimitFlat, upperLimitFlat);

    cout << "interval " << m_lowerLimitOnParent << " " << m_upperLimitOnParent << " flat on " << lowerLimitFlat << " " << upperLimitFlat << " ? "  << isFlat(lowerLimitFlat, upperLimitFlat) << endl;

    uint lowerLimitNewWindow, upperLimitNewWindow;

    if (isFlat(lowerLimitFlat, upperLimitFlat))
    {
        vector<uvec2> roughAreas;
        vector<WLMCWindow::OVERLAPTYPES> overlaps;
        findComplementaryRoughAreas(lowerLimitFlat, upperLimitFlat, roughAreas, overlaps);

        vector<uvec2> _roughAreas(roughAreas.size());
        for (uint i = 0; i < roughAreas.size(); ++i)
        {
            uvec2 roughArea = roughAreas.at(i);
            WLMCWindow::OVERLAPTYPES overlapType = overlaps.at(i);

            getSubWindowLimits(overlapType, roughArea(0), roughArea(1), lowerLimitNewWindow, upperLimitNewWindow);

            _roughAreas.at(i) = {lowerLimitNewWindow, upperLimitNewWindow};
        }

        tmp_output(lowerLimitFlat, upperLimitFlat, _roughAreas);


        for (uint i = 0; i < roughAreas.size(); ++i)
        {
            uvec2 roughArea = roughAreas.at(i);
            WLMCWindow::OVERLAPTYPES overlapType = overlaps.at(i);

            getSubWindowLimits(overlapType, roughArea(0), roughArea(1), lowerLimitNewWindow, upperLimitNewWindow);

            cout << "############### SPAWNED SUB WINDOW ###################### " <<  lowerLimitNewWindow << " " << upperLimitNewWindow << endl;

            m_subWindows.push_back(new WLMCWindow(m_system, m_DOS, m_energies, lowerLimitNewWindow, upperLimitNewWindow, overlapType));
        }

        return getMeanFlatness(lowerLimitFlat, upperLimitFlat);

    }
    else
    {
        tmp_output(0, 0, {{0, m_nbins}});
        return 0;
    }

}


double WLMCWindow::getMeanFlatness(const uint lowerLimit, const uint upperLimit) const
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

void WLMCWindow::findComplementaryRoughAreas(const uint lowerLimitFlat, const uint upperLimitFlat, vector<uvec2> &roughAreas, vector<WLMCWindow::OVERLAPTYPES> &overlaps) const
{

    if (lowerLimitFlat != 0)
    {
        roughAreas.push_back({0, lowerLimitFlat});
        overlaps.push_back(WLMCWindow::OVERLAPTYPES::UPPER);
    }

    if (upperLimitFlat != m_nbins)
    {
        roughAreas.push_back({upperLimitFlat, m_nbins});
        overlaps.push_back(WLMCWindow::OVERLAPTYPES::LOWER);
    }

}

void WLMCWindow::findFlatArea(uint &lowerLimit, uint &upperLimit) const
{
    scanForFlattestArea(lowerLimit, upperLimit);

    expandFlattestArea(lowerLimit, upperLimit);
}

void WLMCWindow::scanForFlattestArea(uint &lowerLimit, uint &upperLimit) const
{
    //performs maximization of flatness for windows of size minWindowSize for increments
    //of size windowIncrementSize. Interval of maximum is stored in lowerLimit and upperLimit.

    uint lowerLimitScan = 0;
    uint upperLimitScan = m_system->minWindowSizeFlat(m_nbins);

    double Fn;

    double maxFn = -1;

    while (upperLimitScan < m_nbins)
    {
        Fn = estimateFlatness(lowerLimitScan, upperLimitScan);

        if (Fn > maxFn)
        {
            maxFn = Fn;

            lowerLimit = lowerLimitScan;
            upperLimit = upperLimitScan;

        }

        lowerLimitScan += m_system->windowIncrementSize();
        upperLimitScan += m_system->windowIncrementSize();
    }

}

void WLMCWindow::expandFlattestArea(uint &lowerLimit, uint &upperLimit) const
{
    //performs maximization of s = u - l under criteria that [u, l] is flat. originates at [lowerLimit, upperLimit]
    //resulting interval where the max of s is found is stored in lowerLimit and upperLimit.

    uint upperLimitStart = upperLimit;

    uint expandingUpperLimit;
    int expandingLowerLimit = lowerLimit;

    uint span;
    uint maxSpan = upperLimit - lowerLimit;
    while (expandingLowerLimit > 0)
    {
        expandingUpperLimit = upperLimitStart;

        while (expandingUpperLimit < m_nbins)
        {
            if (isFlat(expandingLowerLimit, expandingUpperLimit))
            {
                span = expandingUpperLimit - expandingLowerLimit;

                if (span > maxSpan)
                {
                    maxSpan = span;

                    lowerLimit = expandingLowerLimit;
                    upperLimit = expandingUpperLimit;

                }
            }

            expandingUpperLimit += m_system->windowIncrementSize();
        }

        expandingLowerLimit -= m_system->windowIncrementSize();
    }
}

void WLMCWindow::getSubWindowLimits(OVERLAPTYPES overlapType, const uint lowerLimitRough, const uint upperLimitRough, uint &lowerLimit, uint &upperLimit) const
{
    //If the desired area is lower than the minimum window, we need to increase it.
    if (upperLimitRough - lowerLimitRough < m_system->minWindowSizeRough())
    {
        uint halfWindow = m_system->minWindowSizeRough()/2;

        //if expansion on lower side is blocked.
        if (lowerLimitRough < halfWindow)
        {
            lowerLimit = 0;
            upperLimit = m_system->minWindowSizeRough();
        }

        //and the same goes for the upper limit
        else if (upperLimitRough > m_nbins - halfWindow)
        {
            upperLimit = m_nbins;
            lowerLimit = m_nbins - m_system->minWindowSizeRough();
        }

        //if it is not blocked, we simply expand.
        else
        {

            //to account for integer division.
            uint extra = 2*halfWindow - m_system->minWindowSizeRough();
            uint centerPoint = (upperLimitRough + lowerLimitRough)/2;

            upperLimit = centerPoint + halfWindow + extra;
            lowerLimit = centerPoint - halfWindow;

        }

    }

    else
    {
        lowerLimit = lowerLimitRough;
        upperLimit = upperLimitRough;
    }

    //We are not guaranteed a window of size at least minWindow.
    //Now we add the overlap.

    if (overlapType == WLMCWindow::OVERLAPTYPES::LOWER)
    {
        if (lowerLimit < m_system->overlap())
        {
            cout << "ERROR IN LOWER OVERLAP " << lowerLimitRough << " " << upperLimitRough << " " << lowerLimit << endl;
            exit(1);
        }

        lowerLimit -= m_system->overlap();
    }
    else
    {
        if (upperLimit > m_nbins - m_system->overlap())
        {
            cout << "ERROR IN UPPER OVERLAP " << lowerLimitRough << " " << upperLimitRough << " " << upperLimit << endl;
            exit(1);
        }

        upperLimit += m_system->overlap();
    }

}

void WLMCWindow::registerVisit(const uint bin)
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

uint WLMCWindow::getBin(double value) const
{
    uint bin = m_nbins*(value - m_minValue)/m_valueSpan;

    if (bin == m_nbins)
    {
        return m_nbins - 1;
    }

    return bin;
}

void WLMCWindow::reset()
{
    for (WLMCWindow *subWindow : m_subWindows)
    {
        delete subWindow;
    }

    m_subWindows.clear();

    m_visitCounts.fill(m_unsetCount);
}

bool WLMCWindow::isFlat(const uint lowerLimit, const uint upperLimit) const
{
    return estimateFlatness(lowerLimit, upperLimit) >= m_system->flatnessCriterion();
}

void WLMCWindow::mergeWith(WLMCWindow *other, double meanVisitAtFlatArea)
{

    double shiftOldNew = 0;
    double shiftVisit = 0;

    vec oldDOS;
    vec newDOS;

    span overlap;
    span subwindowOverlap;

    if (other->overlapType() == WLMCWindow::OVERLAPTYPES::LOWER)
    {
        overlap = span(other->lowerLimitOnParent(), other->lowerLimitOnParent() + m_system->overlap() - 1);
        subwindowOverlap = span(0, m_system->overlap() - 1);
    }

    else
    {
        overlap = span(other->upperLimitOnParent() - m_system->overlap(), other->upperLimitOnParent() - 1);
        subwindowOverlap = span(other->nbins() - m_system->overlap(), other->nbins() - 1);
    }

    oldDOS = m_DOS(overlap);
    newDOS = other->DOS()(subwindowOverlap);

    shiftOldNew = mean(oldDOS/newDOS);

    m_DOS(overlap) = (oldDOS + newDOS*shiftOldNew)/2;

    cout << "merged " << other->lowerLimitOnParent() << " " << other->upperLimitOnParent() << " with shift on " << shiftOldNew << endl;


    shiftVisit = meanVisitAtFlatArea/other->getMeanFlatness();
    cout << shiftVisit << " " << meanVisitAtFlatArea << " " << other->getMeanFlatness() << endl;
    cout << "shift visit = " << shiftVisit << endl;

    for (uint i = 0; i < other->nbins(); ++i)
    {
        if (!other->isUnsetCount(i) || isUnsetCount(other->lowerLimitOnParent() + i))
        {
            m_visitCounts(other->lowerLimitOnParent() + i) = other->visitCounts(i)*shiftVisit;
        }
    }

}

uint WLMCWindow::topIncrement(const uint upperLimit) const
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

uint WLMCWindow::bottomIncrement(const uint lowerLimit) const
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

void WLMCWindow::tmp_output(const uint lowerLimitFlat, const uint upperLimitFlat, const vector<uvec2> &roughAreas) const
{

    ofstream fFlat;
    fFlat.open(kMC::KMCSolver::instance()->filePath() + "/flatness.txt");

    fFlat << lowerLimitFlat << " " << upperLimitFlat << endl;

    fFlat.close();

    ofstream fRough;
    fRough.open(kMC::KMCSolver::instance()->filePath() + "/roughness.txt");
    for (const uvec2 & roughArea : roughAreas)
    {
        fRough << roughArea(0) << " " << roughArea(1) << endl;
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


void WLMCWindow::adjustLimits(uint &start, uint &end)
{
    double thresh = 1E-16;

    start = 0;
    while (m_DOS(start) <= thresh)
    {
        m_DOS(start) = 0;
        start++;
    }

    end = m_nbins - 1;
    while (m_DOS(end) <= thresh)
    {
        m_DOS(end) = 0;
        end--;
    }

}
