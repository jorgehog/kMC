#include "window.h"
#include "system.h"
#include "windowparams.h"

#include <kMC>
#include <BADAss/badass.h>

using namespace WLMC;

Window::Window(System *system,
               const vec &parentDOS,
               const vec &parentEnergies,
               const vector<bool> &parentBinDeflation,
               const uint lowerLimitOnParent,
               const uint upperLimitOnParent,
               OverlapTypes overlapType,
               bool allowSubwindowing) :
    m_system(system),
    m_allowsSubwindowing(allowSubwindowing),
    m_overlapType(overlapType),
    m_lowerLimitOnParent(lowerLimitOnParent),
    m_upperLimitOnParent(upperLimitOnParent),
    m_nbins(upperLimitOnParent - lowerLimitOnParent),
    m_minValue(parentEnergies(lowerLimitOnParent)),
    m_maxValue(parentEnergies(upperLimitOnParent - 1)),
    m_valueSpan(m_maxValue - m_minValue),
    m_DOS(parentDOS(span(lowerLimitOnParent, upperLimitOnParent - 1))),
    m_energies(parentEnergies(span(lowerLimitOnParent, upperLimitOnParent - 1))),
    m_visitCounts(uvec(m_nbins)),
    m_gradientSampleCounter(0),
    m_spanSum(0),
    m_centerSum(0),
    m_deflatedBins(vector<bool>(m_nbins))
{
    for (uint i = 0; i < m_nbins; ++i)
    {
        m_visitCounts(i) = m_unsetCount;
        m_deflatedBins.at(i) = parentBinDeflation.at(m_lowerLimitOnParent + i);
    }


    BADAss(upperLimitOnParent, >, lowerLimitOnParent, "illegal window.");

}

Window::Window(System *system,
               const uint nBins,
               const double minValue,
               const double maxValue,
               bool allowSubwindowing) :
    Window(system,
           ones(nBins),
           linspace(minValue, maxValue, nBins),
           vector<bool>(nBins, false),
           0,
           nBins,
           Window::OverlapTypes::None,
           allowSubwindowing)
{

}

Window::Window(const Window &parentWindow, const WindowParams &windowParams) :
    Window(parentWindow.system(),
           parentWindow.DOS(),
           parentWindow.energies(),
           parentWindow.deflatedBins(),
           windowParams.m_lowerLimit,
           windowParams.m_upperLimit,
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
    m_system->loadConfigurationForWindow(this);

    BADAssBool(isLegal(m_system->getTotalValue()),
               "Loaded configuration is outside the windowed values.",
               badass::quickie([this] ()
    {
        cout << m_minValue << " " << m_maxValue << endl;
        cout << m_system->getTotalValue() << endl;
    }
    )
               );
}

void Window::calculateWindow()
{
    BADAssBool(isLegal(m_system->getTotalValue()), "Initigal configuration is illegal.", badass::quickie([&] ()
    {
        cout << m_minValue << " " << m_maxValue << " " << m_system->getTotalValue() << endl;
    }));

    cout << "sampling on " << m_lowerLimitOnParent << " " << m_upperLimitOnParent << " f = " << m_system->f() << endl;

    //Need a method for saying that you are flat if flat on the area below overlap. Same goes for continuity.
    while (m_subWindows.empty() && !isFlatOnParent())
    {

        m_system->sampleWindow(this);

        normaliseDOS();

        findSubWindows();

        tmp_output();

        //tmp
        if (!m_subWindows.empty())
        {
            for (Window *w: m_subWindows)
            {
                cout << "subWindow found at " << w->lowerLimitOnParent() << " " << w->upperLimitOnParent() << endl;
            }
        }
        //

        for (Window *subWindow : m_subWindows)
        {
            subWindow->loadInitialConfig();
            subWindow->calculateWindow();
            mergeWith(subWindow);
        }

        tmp_output();
        cout << estimateFlatness() << " " << isFlat() << " " << isFlatOnParent() << endl;
    }

    cout << "Window done: " << m_lowerLimitOnParent << " " << m_upperLimitOnParent << endl;

}

double Window::estimateFlatness(const uint lowerLimit, const uint upperLimit) const
{

    BADAss(upperLimit - lowerLimit, >=, m_system->minWindowSize(),
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
        cout << "continuing: subwindowing illegal." << endl;
        return;
    }

    else if (!findFlatArea())
    {
        cout << "continuing: no flat area found." << endl;
        return;
    }

    else if (!flatProfileIsContinousOnParent())
    {
        cout << "continuing: flat area not contonous on parent." << endl;
        return;
    }

    else if (!flatspanGradientConverged())
    {
        cout << "continuing: flat area has not converged." << endl;
        return;
    }

    else
    {

        cout << "found flat area " << m_flatAreaLower << " " << m_flatAreaUpper << " " << estimateFlatness(m_flatAreaLower, m_flatAreaUpper) << endl;
        tmp_output();

        vector<WindowParams> roughWindowParams;
        findComplementaryRoughAreas(roughWindowParams);

        for (WindowParams &windowParam : roughWindowParams)
        {
            getSubWindowLimits(windowParam);

            m_subWindows.push_back(new Window(*this, windowParam));

        }

    }

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
    if (!overlapsAtBottom())
    {

        if (m_flatAreaLower != 0)
        {
            roughWindowParams.push_back(WindowParams(0, m_flatAreaLower, Window::OverlapTypes::Upper));
        }

    }

    if (!overlapsAtTop())
    {

        if (m_flatAreaUpper != m_nbins)
        {
            roughWindowParams.push_back(WindowParams(m_flatAreaUpper, m_nbins, Window::OverlapTypes::Lower));
        }

    }

}

bool Window::findFlatArea()
{
    if (scanForFlattestArea())
    {
        expandFlattestArea();

        uint center = (m_flatAreaLower + m_flatAreaUpper)/2;
        uint span = m_flatAreaUpper - m_flatAreaLower;

        m_centerSum += center;
        m_spanSum += span;

        m_gradientSampleCounter++;

        for (uint i = 1; i < 4; ++i)
        {
            m_centerSums(i - 1) = m_centerSums(i);
            m_spanSums(i - 1) = m_spanSums(i);
        }

        m_centerSums(3) = m_centerSum/m_gradientSampleCounter;
        m_spanSums(3) = m_spanSum/m_gradientSampleCounter;

        vec weightsDouble = {-1, 4, -5, 2};
        vec weightsSingle = {0, 1.5, -3.0, 1.5};

        m_centerGradient = 0;
        m_spanGradient = 0;
        m_spanLaplace = 0;

        for (uint i = 0; i < 4; ++i)
        {
            m_centerGradient += weightsSingle(i)*m_centerSums(i);
            m_spanGradient += weightsSingle(i)*m_spanSums(i);
            m_spanLaplace += weightsDouble(i)*m_spanSums(i);
        }

        return true;
    }

    else
    {
        m_flatAreaLower = m_flatAreaUpper = 0;
        return false;
    }

}

bool Window::scanForFlattestArea()
{
    //performs maximization of flatness for windows of size minWindowSize for increments
    //of size windowIncrementSize. Interval of maximum is stored in lowerLimit and upperLimit.

    uint lowerLimitScan = 0;
    uint upperLimitScan = m_system->minWindowSize();

    double Fn;

    constexpr double unsetFn = -1;

    double maxFn = unsetFn;

    const double maximalSparsity = 0.5;

    while (upperLimitScan < m_nbins)
    {
        Fn = estimateFlatness(lowerLimitScan, upperLimitScan);

        if (Fn > maxFn)
        {
            if (sparsity(lowerLimitScan, upperLimitScan) < maximalSparsity)
            {

                maxFn = Fn;

                m_flatAreaLower = lowerLimitScan;
                m_flatAreaUpper = upperLimitScan;

            }

        }

        lowerLimitScan++;
        upperLimitScan++;
    }

    return maxFn >= m_system->flatnessCriterion();

}

double Window::sparsity(const uint lowerLimit, const uint upperLimit) const
{
    double sparsity = 0;

    for (uint bin = lowerLimit; bin < upperLimit; ++bin)
    {
        if (isDeflatedBin(bin) || isUnsetCount(bin))
        {
            sparsity++;
        }
    }

    return sparsity/(upperLimit - lowerLimit);
}

void Window::expandFlattestArea()
{
    //performs maximization of s = u - l under criteria that [u, l] is flat. originates at [lowerLimit, upperLimit]
    //resulting interval where the max of s is found is stored in lowerLimit and upperLimit.

    int expandingLowerLimit = m_flatAreaLower;
    uint upperLimitStart = m_flatAreaUpper;

    uint maxSpan = m_flatAreaUpper - m_flatAreaLower;

    while (expandingLowerLimit >= 0)
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

            expandingUpperLimit++;
        }

        expandingLowerLimit--;
    }

}

bool Window::flatProfileIsContinousOnParent() const
{
    if (overlapsAtBottom())
    {
        cout << m_flatAreaLower << " < " << m_system->overlap() << " ? " << endl;
        return m_flatAreaLower <= m_system->overlap();
    }

    else if (overlapsAtTop())
    {
        cout << m_flatAreaUpper << " > " << m_nbins - m_system->overlap() << " ? " << endl;
        return m_flatAreaUpper >= m_nbins - m_system->overlap();
    }

    else
    {
        return true;
    }
}

void Window::getSubWindowLimits(WindowParams &windowParams) const
{

    uint span = windowParams.m_upperLimit - windowParams.m_lowerLimit;

    //Expand upper or lower limits that are too small. This disables the new window from having children.
    if (span < m_system->minWindowSize())
    {
        cout << "span " << span << " too low for window size " << m_system->minWindowSize() << endl;

        windowParams.m_allowSubwindowing = false;

        if (windowParams.m_overlapType == Window::OverlapTypes::Lower)
        {
            BADAss(windowParams.m_upperLimit, >=, m_system->minWindowSize());

            windowParams.m_lowerLimit = windowParams.m_upperLimit - m_system->minWindowSize();
        }

        else if (windowParams.m_overlapType == Window::OverlapTypes::Upper)
        {
            windowParams.m_upperLimit = windowParams.m_lowerLimit + m_system->minWindowSize();

            BADAss(windowParams.m_upperLimit, <=, m_nbins);
        }

    }


    //Add the overlap
    if (windowParams.m_overlapType == Window::OverlapTypes::Lower)
    {
        BADAss(windowParams.m_lowerLimit, >=, m_system->overlap());

        windowParams.m_lowerLimit -= m_system->overlap();
    }

    else if (windowParams.m_overlapType == Window::OverlapTypes::Upper)
    {
        BADAss(m_nbins - windowParams.m_upperLimit, >=, m_system->overlap());

        windowParams.m_upperLimit += m_system->overlap();
    }

}

void Window::registerVisit(const uint bin)
{
    BADAssBool(!isDeflatedBin(bin), "deflated bins should not be visited.");

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


    m_gradientSampleCounter = 0;

    m_spanSum = 0;
    m_centerSum = 0;

    m_subWindows.clear();

    m_visitCounts.fill(m_unsetCount);
}

void Window::deflateDOS()
{
    double eps = 1E-16;
    for (uint bin = 0; bin < m_nbins; ++bin)
    {

        if (isDeflatedBin(bin))
        {
            continue;
        }

        else if (m_DOS(bin) < eps)
        {
            m_visitCounts(bin) = m_unsetCount;
            m_DOS(bin) = 0;
            m_deflatedBins.at(bin) = true;
        }
    }
}

void Window::resetDOS()
{
    m_DOS.ones();
}

bool Window::allowsSubwindowing() const
{
    return m_allowsSubwindowing;
}

bool Window::flatspanGradientConverged() const
{
    if (m_gradientSampleCounter < m_centerSums.n_elem)
    {
        return false;
    }

    const double lim = 0.1;

    cout << "CG " << m_centerGradient << endl;
    cout << "SG " << m_spanGradient << endl;
    cout << "SL " << m_spanLaplace << endl;

    return m_spanLaplace <= 0  && fabs(m_spanGradient) < lim && fabs(m_centerGradient) < lim;
}

bool Window::isFlat(const uint lowerLimit, const uint upperLimit) const
{
    return estimateFlatness(lowerLimit, upperLimit) >= m_system->flatnessCriterion();
}

bool Window::isFlatOnParent() const
{
    if (overlapsAtTop())
    {
        return isFlat(0, m_nbins - m_system->overlap());
    }
    else if (overlapsAtBottom())
    {
        return isFlat(m_system->overlap(), m_nbins);
    }
    else
    {
        return isFlat();
    }
}

void Window::mergeWith(Window *other)
{
    cout << m_lowerLimitOnParent << " " << m_upperLimitOnParent << " ############### MERGED SUB WINDOW ###################### " <<  other->lowerLimitOnParent() << " " << other->upperLimitOnParent() << endl;

    span spanOnParent, _span;
    uint overlapPointOnParent, overlapPoint;

    other->normaliseDOS();
    normaliseDOS();


    overlapPointOnParent = getOverlapPoint(other);
    if (overlapPointOnParent == m_unsetCount)
    {
        return;
    }

    overlapPoint = overlapPointOnParent - other->lowerLimitOnParent();

    if (other->overlapsAtBottom())
    {
        _span = span(overlapPoint, other->nbins() - 1);
        spanOnParent = span(overlapPointOnParent, other->upperLimitOnParent() - 1);
    }

    else
    {
        _span = span(0, overlapPoint - 1);
        spanOnParent = span(other->lowerLimitOnParent(), overlapPointOnParent - 1);
    }

    double shiftDOS = m_DOS(overlapPointOnParent)/other->DOS(overlapPoint);

    m_DOS(spanOnParent) = shiftDOS*other->DOS()(_span);

    normaliseDOS();

}

uint Window::getOverlapPoint(const Window *other)
{
    uint init;
    int bin;

    if (other->overlapsAtBottom())
    {
        init = (m_flatAreaUpper + other->lowerLimitOnParent())/2;
    }

    else
    {
        init = (m_flatAreaLower + other->upperLimitOnParent())/2;
    }

    bin = init;

    while (isDeflatedBin(bin) || other->isDeflatedBin(bin - other->lowerLimitOnParent()))
    {

        if (other->overlapsAtBottom())
        {

            bin++;

            if (bin >= signed(other->upperLimitOnParent()))
            {
                m_DOS(span(init, other->upperLimitOnParent() - 1)).zeros();
                return m_unsetCount;
            }

        }

        else
        {
            bin--;

            if (bin < signed(other->lowerLimitOnParent()))
            {
                m_DOS(span(other->lowerLimitOnParent(), init - 1)).zeros();
                return m_unsetCount;
            }

        }

    }

    return bin;

}

void Window::tmp_output() const
{

    ofstream fFlat;
    fFlat.open(kMC::KMCSolver::instance()->filePath() + "/flatness.txt");

    fFlat << m_flatAreaLower << " " << m_flatAreaUpper << endl;

    fFlat.close();

    ofstream fRough;
    fRough.open(kMC::KMCSolver::instance()->filePath() + "/roughness.txt");
    for (const Window *subWindow: m_subWindows)
    {
        fRough << subWindow->lowerLimitOnParent() << " " << subWindow->upperLimitOnParent() << endl;
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
