#include "wlmcwindow.h"
#include "wlmcsystem.h"

using namespace WLMC;

WLMCWindow::WLMCWindow(WLMCSystem *system, const vec &DOS,
                       const uint lowerLimit,
                       const uint upperLimit,
                       const double minValue,
                       const double maxValue) :
    m_system(system),
    m_lowerLimit(lowerLimit),
    m_upperLimit(upperLimit),
    m_nbins(upperLimit - lowerLimit),
    m_minValue(minValue),
    m_maxValue(maxValue),
    m_valueSpan(maxValue - minValue),
    m_DOS(DOS),
    m_visitCounts(uvec(m_nbins))
{
    m_visitCounts.fill(m_unsetCount);
}

WLMCWindow::WLMCWindow(WLMCSystem *system, const uint nBins,
                       const double minValue,
                       const double maxValue) :
    WLMCWindow(system, ones(nBins), 0, nBins, minValue, maxValue)
{

}

WLMCWindow::~WLMCWindow()
{
    m_DOS.clear();
    m_visitCounts.clear();
    m_system = NULL;
}

void WLMCWindow::calculateWindow(kMC::KMCSolver *solver)
{
    m_system->sampleWindow(this);

    m_DOS = normalise(m_DOS);

    vector<uvec2> flatAreas;
    vector<uvec2> roughAreas;

    findFlatAreas(flatAreas, m_lowerLimit, m_upperLimit);
    findComplementaryRoughAreas(flatAreas, roughAreas, m_lowerLimit, m_upperLimit);

    //todo: definine flatness minimum window size by the size of the parent window. i.e. must be flat on atleast 50%.
    //todo: debug singleWindow version.
    //todo: initialize subwindows.
    //todo: merge

    tmp_output(solver, flatAreas, roughAreas);
    
}

double WLMCWindow::estimateFlatness(const uvec &visitCounts) const
{
    uint min = std::numeric_limits<uint>::max();

    double mean = 0;

    uint nSetCounts = 0;

    for (const uint &vc : visitCounts)
    {
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

void WLMCWindow::findFlatAreas(vector<uvec2> &flatAreas, const uint lowerLimit, const uint upperLimit) const
{
    uint upperLimitFlat, lowerLimitFlat, origin;

    origin = findFlattestOrigin(lowerLimit, upperLimit);

    findFlatArea(upperLimitFlat, lowerLimitFlat, origin);

    if (isFlat(m_visitCounts(span(lowerLimitFlat, upperLimitFlat - 1))))
    {
        flatAreas.push_back(uvec2({lowerLimitFlat, upperLimitFlat}));

        vector<uvec2> complementaryRoughAreas;
        findComplementaryRoughAreas(flatAreas, complementaryRoughAreas, lowerLimit, upperLimit);

        for (const uvec2 &complementaryRoughArea : complementaryRoughAreas)
        {
            findFlatAreas(flatAreas, complementaryRoughArea(0), complementaryRoughArea(1));
        }

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

void WLMCWindow::findComplementaryRoughAreas(const vector<uvec2> &flatAreas, vector<uvec2> &roughAreas, const uint lowerLimit, const uint upperLimit) const
{
    if (flatAreas.empty())
    {
        roughAreas.push_back({lowerLimit, upperLimit});
        return;
    }

    uint prev = lowerLimit;

    for (const uvec2 &flatArea : flatAreas)
    {
        if (prev != flatArea(0))
        {
            roughAreas.push_back({prev, flatArea(0)});
        }

        prev = flatArea(1);
    }

    uint endLastFlat = flatAreas.back()(1);

    if (endLastFlat != upperLimit)
    {
        roughAreas.push_back({endLastFlat, upperLimit});
    }

}

void WLMCWindow::findFlatArea(uint &upperLimit, uint &lowerLimit, const uint origin) const
{

    if (origin > m_nbins - m_system->minWindowSize()/2)
    {
        upperLimit = m_nbins;
        lowerLimit = origin - m_system->minWindowSize();
    }

    else if (origin < m_system->minWindowSize()/2)
    {
        lowerLimit = 0;
        upperLimit = m_system->minWindowSize();
    }

    else
    {
        upperLimit = origin + m_system->minWindowSize()/2;
        lowerLimit = origin - m_system->minWindowSize()/2;
    }

    bool upperSet = false;
    bool lowerSet = false;

    //right propagation
    while (!(upperSet && lowerSet))
    {

        if (!upperSet)
        {

            if (upperLimit == m_nbins)
            {
                upperSet = true;
            }

            if (!isFlat(m_visitCounts(span(lowerLimit, topIncrement(upperLimit) - 1))))
            {
                upperSet = true;
            }
            else
            {
                upperLimit = topIncrement(upperLimit);
            }

        }

        if (!lowerSet)
        {

            if (lowerLimit == 0)
            {
                lowerSet = true;
            }

            if (!isFlat(m_visitCounts(span(bottomIncrement(lowerLimit), upperLimit - 1))))
            {
                lowerSet = true;
            }
            else
            {
                lowerLimit = bottomIncrement(lowerLimit);
            }

        }

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
    m_visitCounts.fill(m_unsetCount);
}

bool WLMCWindow::isFlat(const uvec &visitCounts) const
{
    return estimateFlatness(visitCounts) >= m_system->flatnessCriterion();
}

void WLMCWindow::mergeWith(WLMCWindow *other)
{

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

void WLMCWindow::tmp_output(kMC::KMCSolver *solver, const vector<uvec2> &flatAreas, const vector<uvec2> &roughAreas) const
{

    ofstream fFlat;
    fFlat.open(solver->filePath() + "/flatness.txt");
    for (const uvec2 & flatArea : flatAreas)
    {
        cout << estimateFlatness(m_visitCounts(span(flatArea(0), flatArea(1) - 1))) << " flat on " << flatArea.t();

        fFlat << flatArea(0) << " " << flatArea(1) << endl;
    }
    fFlat.close();

    ofstream fRough;
    fRough.open(solver->filePath() + "/roughness.txt");
    for (const uvec2 & roughArea : roughAreas)
    {
        cout << estimateFlatness(m_visitCounts(span(roughArea(0), roughArea(1) - 1))) << " rough on " << roughArea.t();

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


