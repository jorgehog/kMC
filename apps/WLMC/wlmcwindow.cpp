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
    findFlatAreas(flatAreas);

    //todo: smarter way to generate flat areas. Start where it is flattest.
    //todo: debug singleWindow version.
    //todo: locate flat areas and initialize subwindows.

    tmp_output(solver, flatAreas);
    
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

void WLMCWindow::findFlatAreas(vector<uvec2> &flatAreas) const
{
    uint upperLimit, lowerLimit, origin;

    origin = 0;
    while (origin < m_nbins)
    {
        findFlatArea(upperLimit, lowerLimit, origin);

        if (!isFlat(m_visitCounts(span(lowerLimit, upperLimit - 1))))
        {
            origin += m_system->windowIncrementSize();
        }

        else
        {
            flatAreas.push_back(uvec2({lowerLimit, upperLimit}));
            origin = upperLimit;
        }

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

void WLMCWindow::tmp_output(kMC::KMCSolver *solver, const vector<uvec2> &flatAreas) const
{

    ofstream f;
    f.open(solver->filePath() + "/flatness.txt");
    for (const uvec2 & flatArea : flatAreas)
    {
        cout << estimateFlatness(m_visitCounts(span(flatArea(0), flatArea(1) - 1))) << " flat on " << flatArea.t();

        f << flatArea(0) << " " << flatArea(1) << endl;
    }
    f.close();

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


