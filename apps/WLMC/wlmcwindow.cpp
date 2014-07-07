#include "wlmcwindow.h"
#include "wlmcsystem.h"

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

void WLMCWindow::calculateWindow(kMC::KMCSolver *solver)
{
    cout << "sampling on " << m_lowerLimitOnParent << " " << m_upperLimitOnParent << endl;

    uint lowerLimitFlat, upperLimitFlat;
    vector<uvec2> roughAreas;

    while (!isFlat() || m_subWindows.size() != 0)
    {

        m_system->sampleWindow(this);

        m_DOS = normalise(m_DOS);

        findFlatArea(lowerLimitFlat, upperLimitFlat);
        findComplementaryRoughAreas(lowerLimitFlat, upperLimitFlat, roughAreas);

        tmp_output(solver, lowerLimitFlat, upperLimitFlat, roughAreas);

        roughAreas.clear();

//        findSubWindows(solver);

//        for (WLMCWindow *subWindow : m_subWindows)
//        {
//            subWindow->calculateWindow(solver);
//            mergeWith(subWindow);
//        }

    }

    //todo: new improved robust scanning flat area algorithm.
    //      1. iterate over n windows of size W and make F(n)
    //      2. pick limits which maximizes F(n)
    //      3. find longest interval [n - dn_l, n + dn_u] which is flat
    //todo: debug singleWindow version.
    //todo: initialize subwindows.
    //todo: merge
}

double WLMCWindow::estimateFlatness(const uint lowerLimit, const uint upperLimit) const
{
    if (upperLimit - lowerLimit < m_system->minWindowSize())
    {
        cout << "error. Flatness requested on window lower than minimum limit." << endl;
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

void WLMCWindow::findSubWindows(kMC::KMCSolver *solver)
{
    uint lowerLimitNewWindow, upperLimitNewWindow, upperLimitFlat, lowerLimitFlat;
    OVERLAPTYPES overlapType;

    findFlatArea(lowerLimitFlat, upperLimitFlat);

    cout << "interval " << m_lowerLimitOnParent << " " << m_upperLimitOnParent << " flat on " << lowerLimitFlat << " " << upperLimitFlat << " ? "  << isFlat(lowerLimitFlat, upperLimitFlat) << endl;

    if (isFlat(lowerLimitFlat, upperLimitFlat))
    {
        vector<uvec2> roughAreas;
        findComplementaryRoughAreas(lowerLimitFlat, upperLimitFlat, roughAreas);

        tmp_output(solver, lowerLimitFlat, upperLimitFlat, roughAreas);

        for (const uvec2 roughArea : roughAreas)
        {
            if (roughArea(1) == m_lowerLimitOnParent)
            {
                overlapType = OVERLAPTYPES::UPPER;
            }
            else if (roughArea(0) == m_upperLimitOnParent)
            {
                overlapType = OVERLAPTYPES::LOWER;
            }
            else
            {
                cout << "error" << endl;
                exit(1);
            }

            getSubWindowLimits(overlapType, roughArea(0), roughArea(1), lowerLimitNewWindow, upperLimitNewWindow);
            m_subWindows.push_back(new WLMCWindow(m_system, m_DOS, m_energies, lowerLimitNewWindow, upperLimitNewWindow, overlapType));
        }

    }
    else
    {
        tmp_output(solver, 0, 0, {});
    }

}

uint WLMCWindow::findFlattestOrigin(const uint lowerLimit, const uint upperLimit) const
{

    double mean = getMeanFlatness(lowerLimit, upperLimit);

    uint flattest;
    double minError = std::numeric_limits<double>::max();

    for (uint i = lowerLimit; i < upperLimit; ++i)
    {
        uint vc = m_visitCounts(i);

        if (!(vc == m_unsetCount))
        {

            double error = (vc - mean)*(vc - mean);

            if (error < minError)
            {
                minError = error;
                flattest = i;
            }

        }
    }

    return flattest;

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

void WLMCWindow::findComplementaryRoughAreas(const uint lowerLimitFlat, const uint upperLimitFlat, vector<uvec2> &roughAreas) const
{

    if (lowerLimitFlat != 0)
    {
        roughAreas.push_back({0, lowerLimitFlat});
    }

    if (upperLimitFlat != m_nbins)
    {
        roughAreas.push_back({upperLimitFlat, m_nbins});
    }

}

void WLMCWindow::findFlatArea(uint &lowerLimit, uint &upperLimit) const
{

    scanForFlattestArea(lowerLimit, upperLimit);

    expandFlattestArea(lowerLimit, upperLimit);

//    cout << "origin : " << origin << endl;

//    if (origin > m_nbins - m_system->minWindowSize()/2)
//    {
//        upperLimit = m_nbins;
//        lowerLimit = origin - m_system->minWindowSize();
//    }

//    else if (origin < m_system->minWindowSize()/2)
//    {
//        lowerLimit = 0;
//        upperLimit = m_system->minWindowSize();
//    }

//    else
//    {
//        upperLimit = origin + m_system->minWindowSize()/2;
//        lowerLimit = origin - m_system->minWindowSize()/2;
//    }

//    bool upperSet = false;
//    bool lowerSet = false;

//    //right propagation
//    while (!(upperSet && lowerSet))
//    {

//        if (!upperSet)
//        {

//            if (upperLimit == m_nbins)
//            {
//                upperSet = true;
//            }

//            if (!isFlat(m_visitCounts(span(lowerLimit, topIncrement(upperLimit) - 1))))
//            {
//                upperSet = true;
//            }
//            else
//            {
//                upperLimit = topIncrement(upperLimit);
//            }

//        }

//        if (!lowerSet)
//        {

//            if (lowerLimit == 0)
//            {
//                lowerSet = true;
//            }

//            if (!isFlat(m_visitCounts(span(bottomIncrement(lowerLimit), upperLimit - 1))))
//            {
//                lowerSet = true;
//            }
//            else
//            {
//                lowerLimit = bottomIncrement(lowerLimit);
//            }

//        }

    //    }
}

void WLMCWindow::scanForFlattestArea(uint &lowerLimit, uint &upperLimit) const
{
    //performs maximization of flatness for windows of size minWindowSize for increments
    //of size windowIncrementSize. Interval of maximum is stored in lowerLimit and upperLimit.

    uint lowerLimitScan = 0;
    uint upperLimitScan = m_system->minWindowSize();

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
    uint maxSpan = 0;
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
    if (upperLimitRough - lowerLimitRough < m_system->minWindowSize())
    {

        uint halfWindow = m_system->minWindowSize()/2;

        //if expansion on lower side is blocked.
        if (lowerLimitRough < m_lowerLimitOnParent + halfWindow)
        {
            lowerLimit = m_lowerLimitOnParent;
            upperLimit = m_lowerLimitOnParent + m_system->minWindowSize();
        }

        //and the same goes for the upper limit
        else if (upperLimitRough > m_upperLimitOnParent - halfWindow)
        {
            upperLimit = m_upperLimitOnParent;
            lowerLimit = m_upperLimitOnParent - m_system->minWindowSize();
        }

        //if it is not blocked, we simply expand.
        else
        {

            //to account for integer division.
            uint extra = 2*halfWindow - m_system->minWindowSize();
            uint centerPoint = (upperLimitRough + lowerLimitRough)/2;

            upperLimit = centerPoint + halfWindow + extra;
            lowerLimit = centerPoint - halfWindow;

        }

    }

    //We are not guaranteed a window of size at least minWindow.
    //Now we add the overlap.

    if (overlapType == WLMCWindow::OVERLAPTYPES::LOWER)
    {
        if (lowerLimit < m_system->overlap())
        {
            cout << "ERROR IN OVERLAP" << lowerLimit << endl;
            exit(1);
        }

        lowerLimit -= m_system->overlap();
    }
    else
    {
        if (upperLimit > m_nbins - m_system->overlap())
        {
            cout << "ERROR IN OVERLAP" << upperLimit << endl;
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

void WLMCWindow::mergeWith(WLMCWindow *other)
{
    m_DOS(span(other->lowerLimit(), other->upperLimit() - 1)) = other->DOS();

    m_visitCounts(span(other->lowerLimit(), other->upperLimit() - 1)) = other->visitCounts();
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

void WLMCWindow::tmp_output(kMC::KMCSolver *solver, const uint lowerLimitFlat, const uint upperLimitFlat, const vector<uvec2> &roughAreas) const
{

    ofstream fFlat;
    fFlat.open(solver->filePath() + "/flatness.txt");

    if (lowerLimitFlat != 0 && upperLimitFlat != 0)
    {
        cout << estimateFlatness(lowerLimitFlat, upperLimitFlat) << " flat on " << lowerLimitFlat << " " << upperLimitFlat << endl;
    }

    fFlat << lowerLimitFlat << " " << upperLimitFlat << endl;

    fFlat.close();

    ofstream fRough;
    fRough.open(solver->filePath() + "/roughness.txt");
    for (const uvec2 & roughArea : roughAreas)
    {
        cout << estimateFlatness(roughArea(0), roughArea(1)) << " rough on " << roughArea.t();

        fRough << roughArea(0) << " " << roughArea(1) << endl;
    }
    fRough.close();


    stringstream file;
    file << solver->filePath() << "/stateDensity" << m_system->nParticles() << ".arma";
    vec E = linspace<vec>(m_minValue, m_maxValue, m_nbins);

    uvec indices = find(m_visitCounts != m_unsetCount);
    uvec vc = m_visitCounts(indices);

    vec dos = m_DOS(indices);
    vec e = E(indices);
    vec idx = conv_to<vec>::from(indices);
    vec vcd = conv_to<vec>::from(vc);

    if (join_rows(join_rows(join_rows(e, dos), vcd), idx).eval().save(file.str()))
    {
        cout << "storing " << file.str() << endl;
    }
}


