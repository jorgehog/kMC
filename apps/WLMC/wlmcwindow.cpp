#include "wlmcwindow.h"
#include "wlmcsystem.h"

using namespace WLMC;

WLMCWindow::WLMCWindow(const vec &DOS,
                       const uint lowerLimit,
                       const uint upperLimit,
                       const double minValue,
                       const double maxValue) :
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

void WLMCWindow::calculateWindow(double f)
{
    m_f = f;


}

double WLMCWindow::estimateFlatness() const
{
    uint min = std::numeric_limits<uint>::max();

    double mean = 0;

    uint nSetCounts = 0;

    for (uint i = m_lowerLimit; i < m_upperLimit; ++i)
    {
        if (!isUnsetCount(i))
        {
            uint count = m_visitCounts(i);

            if (count < min)
            {
                count = min;
            }

            mean += count;

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

    m_DOS(bin) *= m_f;
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

void WLMCWindow::performMove()
{

}

void WLMCWindow::mergeWith(WLMCWindow *other)
{

}
