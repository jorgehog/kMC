#include "system.h"
#include "window.h"
#include "windowparams.h"

#include <kMC>
#include <BADAss/badass.h>

using namespace WLMC;

System::System(const uint nParticles,
               const uint NX,
               const uint NY,
               const uint NZ,
               const uint movesPerSampling,
               const double flatnessCriterion,
               const uint overlap,
               const uint minWindowSizeRough,
               function<double()> URNG) :
    m_nParticles(nParticles),
    m_NX(NX),
    m_NY(NY),
    m_NZ(NZ),
    m_volume(m_NX*m_NY*m_NZ),
    m_freeVolume(m_volume - nParticles),
    m_movesPerSampling(movesPerSampling),
    m_flatnessCriterion(flatnessCriterion),
    m_overlap(overlap),
    m_minWindowSize(minWindowSizeRough),
    m_URNG(URNG)
{

    BADAss(nParticles, >=, 2, "Number of particles.");
    BADAss(nParticles, <, m_volume, "Incorrect number of particles.");

}

void System::execute(const uint nbins, const double adaptive, const double fStart, const double fEnd, function<double(double)> reduceFunction)
{

    double max, min;
    locateGlobalExtremaValues(min, max);

    Window mainWindow(this, nbins, min, max, adaptive);

    //    clipWindow(mainWindow);

    setupPresetWindowConfigurations(mainWindow);

//    Window subWindow(mainWindow, WindowParams(0, 150, Window::OverlapTypes::Upper));

//    subWindow.loadInitialConfig();
//    m_f = fStart;
//    cout << " alal " << getTotalValue() << endl;
//    subWindow.calculateWin|dow();
//    return;

    m_f = fStart;
    while (m_f >= fEnd)
    {
        mainWindow.calculateWindow();

        m_f = reduceFunction(m_f);

        mainWindow.reset();
    }

}

void System::sampleWindow(Window *window)
{
    uint nMoves = 0;
    
    while (nMoves < m_movesPerSampling)
    {
        if (doWLMCMove(window))
        {
            nMoves++;
        }
    }
}

bool System::doWLMCMove(Window *window)
{
    uint particleIndex, xd, yd, zd;

    getRandomParticleAndDestination(particleIndex, xd, yd, zd);

    double oldValue = getTotalValue();
    uint oldBin = window->getBin(oldValue);

    double newValue = oldValue + getValueDifference(particleIndex, xd, yd, zd);

    if (!window->isLegal(newValue))
    {
        return false;
    }

    uint newBin = window->getBin(newValue);

    if (window->isDeflatedBin(oldBin) || window->isDeflatedBin(newBin))
    {
        changePosition(particleIndex, xd, yd, zd);
        return false;
    }


    double oldDOS = window->DOS(oldBin);
    double newDOS = window->DOS(newBin);

    BADAss(newDOS, >, 1E-16);

    bool accepted = true;

    if (oldDOS < newDOS)
    {
        accepted = (m_URNG() < oldDOS/newDOS);
    }

    if (accepted)
    {
        changePosition(particleIndex, xd, yd, zd);

        window->registerVisit(newBin);
    }

    else
    {
        window->registerVisit(oldBin);
    }

    return true;

}

void System::doRandomMove()
{
    uint particleIndex, xd, yd, zd;

    xd = yd = zd = 0;

    getRandomParticleAndDestination(particleIndex, xd, yd, zd);
    changePosition(particleIndex, xd, yd, zd);
}

void System::findDestination(const uint destination, uint &xd, uint &yd, uint &zd) const
{
    uint search = 0;

    xd = yd = zd = 0;

    for (uint x = 0; x < m_NX; ++x)
    {
        for (uint y = 0; y < m_NY; ++y)
        {
            for (uint z = 0; z < m_NZ; ++z)
            {
                if (!isOccupiedLoction(x, y, z))
                {
                    if (search == destination)
                    {
                        xd = x;
                        yd = y;
                        zd = z;

                        return;
                    }

                    search++;
                }
            }
        }
    }

}

void System::locateGlobalExtremaValues(double &min, double &max)
{
    uint nSweeps, sweep, x, y, z;
    double localMax, localMin;
    bool isIn;
    vector<double> allExtrema;

    nSweeps = 1;
    sweep = 0;
    min = std::numeric_limits<double>::max();
    m_presetMinimum.set_size(m_nParticles, 3);

    while (sweep < nSweeps)
    {
        randomizeParticlePositions();

        localMin = getGlobalExtremum(System::extrema::minimum);

        isIn = false;
        for (double extrema : allExtrema)
        {
            if (fabs(localMin - extrema) < 0.001)
            {
                isIn = true;
                break;
            }
        }
        if (!isIn)
        {
            allExtrema.push_back(localMin);
            kMC::KMCSolver::instance()->dumpLAMMPS(allExtrema.size());
        }

        if (localMin < min)
        {
            min = localMin;

            for (uint particleIndex = 0; particleIndex < m_nParticles; ++particleIndex)
            {
                getPosition(particleIndex, x, y, z);

                m_presetMinimum(particleIndex, 0) = x;
                m_presetMinimum(particleIndex, 1) = y;
                m_presetMinimum(particleIndex, 2) = z;
            }
        }

        sweep++;
    }

    nSweeps = 20;
    sweep = 0;
    max = std::numeric_limits<double>::min();
    m_presetMaximum.set_size(m_nParticles, 3);

    while (sweep < nSweeps)
    {
        randomizeParticlePositions();

        localMax = getGlobalExtremum(System::extrema::maximum);

        isIn = false;
        for (double extrema : allExtrema)
        {
            if (fabs(localMax - extrema) < 0.001)
            {
                isIn = true;
                break;
            }
        }
        if (!isIn)
        {
            allExtrema.push_back(localMax);
            kMC::KMCSolver::instance()->dumpLAMMPS(allExtrema.size());
        }

        if (localMax > max)
        {
            max = localMax;

            for (uint particleIndex = 0; particleIndex < m_nParticles; ++particleIndex)
            {
                getPosition(particleIndex, x, y, z);

                m_presetMaximum(particleIndex, 0) = x;
                m_presetMaximum(particleIndex, 1) = y;
                m_presetMaximum(particleIndex, 2) = z;
            }
        }

        sweep++;
    }

    cout << "found " << allExtrema.size() << " extrema." << endl;

}

void System::setupPresetWindowConfigurations(const Window &mainWindow)
{
    uint x, y, z, bin, nSet;
    double value, max, min;

    max = mainWindow.maxValue();
    min = mainWindow.minValue();

    uint n = ceil(mainWindow.nbins()/double(m_minWindowSize + m_overlap));

    m_presetWindowConfigurations.set_size(n, m_nParticles, 3);
    m_presetWindowValues = linspace(min, max, n + 1);

    uvec binSet = zeros<uvec>(n);

    //Set the max and min as extrema configurations
    m_presetWindowConfigurations(span(0), span::all, span::all) = m_presetMinimum;
    m_presetWindowConfigurations(span(n-1), span::all, span::all) = m_presetMaximum;

    binSet(0) = 1;
    binSet(n - 1) = 1;
    nSet = 2;

    while (nSet < n)
    {
        doRandomMove();

        value = getTotalValue();

        if (value < min || value > max)
        {
            continue;
        }

        bin = getPresetBinFromValue(value);

        if (binSet(bin) != 0)
        {
            continue;
        }

        binSet(bin) = 1;

        for (uint particleIndex = 0; particleIndex < m_nParticles; ++particleIndex)
        {
            getPosition(particleIndex, x, y, z);

            m_presetWindowConfigurations(bin, particleIndex, 0) = x;
            m_presetWindowConfigurations(bin, particleIndex, 1) = y;
            m_presetWindowConfigurations(bin, particleIndex, 2) = z;
        }

        cout << "set energy for value " << value << " in bin " << bin << " with energyspan " << m_presetWindowValues(bin) << " - " << m_presetWindowValues(bin + 1) << endl;

        nSet++;
    }
}

void System::loadConfigurationForWindow(const Window *window)
{

    double value = (window->minValue() + window->maxValue())/2;

    uint bin = getPresetBinFromValue(value);

    cout << "value " << value << " belongs in bin " << bin << " between " << m_presetWindowValues(bin) << " - " << m_presetWindowValues(bin + 1) << endl;

    uint xPreset, yPreset, zPreset, x, y, z, xAvailable, yAvailable, zAvailable, particleIndexConfig;

    xAvailable = yAvailable = zAvailable = 0;

    bool isAlreadyOccupied;

    kMC::KMCSolver::instance()->clearParticles();

    for (uint particleIndex = 0; particleIndex < m_nParticles; ++particleIndex)
    {
        kMC::KMCSolver::instance()->forceSpawnParticle(m_presetWindowConfigurations(bin, particleIndex, 0),
                                                       m_presetWindowConfigurations(bin, particleIndex, 1),
                                                       m_presetWindowConfigurations(bin, particleIndex, 2));
    }
    return;


    for (uint particleIndex = 0; particleIndex < m_nParticles; ++particleIndex)
    {
        getPosition(particleIndex, x, y, z);

        isAlreadyOccupied = false;
        for (particleIndexConfig = 0; particleIndexConfig < m_nParticles; ++particleIndexConfig)
        {
            xPreset = m_presetWindowConfigurations(bin, particleIndexConfig, 0);
            yPreset = m_presetWindowConfigurations(bin, particleIndexConfig, 1);
            zPreset = m_presetWindowConfigurations(bin, particleIndexConfig, 2);

            if (xPreset == x && yPreset == y && zPreset == z)
            {
                isAlreadyOccupied = true;
                break;
            }

            else if (!isOccupiedLoction(xPreset, yPreset, zPreset))
            {
                xAvailable = xPreset;
                yAvailable = yPreset;
                zAvailable = zPreset;
            }
        }


        if (!isAlreadyOccupied)
        {
            changePosition(particleIndex, xAvailable, yAvailable, zAvailable);
        }
    }

}

uint System::getPresetBinFromValue(const double value) const
{
    uint bin = 0;
    while (true)
    {
        if (value >= m_presetWindowValues(bin) && value <= m_presetWindowValues(bin + 1))
        {
            return bin;
        }

        bin++;
    }

}

void System::clipWindow(Window &window) const
{

    uint histSamples, upperLimit, lowerLimit;

    if (m_movesPerSampling > 10000)
    {
        histSamples = 10000;
    }
    else
    {
        histSamples = m_movesPerSampling;
    }

    vec hist = window.getHistogram(histSamples);

    double m = arma::max(hist);

    uint clipSize = 10;
    double thresh = 1E-2;

    uint upperLimitClip = window.nbins();
    uint lowerLimitClip = window.nbins() - clipSize;

    while (mean(hist(span(lowerLimitClip, upperLimitClip - 1))) < thresh*m)
    {
        upperLimitClip -= clipSize;
        lowerLimitClip -= clipSize;
    }

    upperLimit = (lowerLimitClip + upperLimitClip)/2;

    lowerLimitClip = 0;
    upperLimitClip = clipSize;

    while (mean(hist(span(lowerLimitClip, upperLimitClip - 1))) < thresh*m)
    {
        upperLimitClip += clipSize;
        lowerLimitClip += clipSize;
    }

    lowerLimit = (lowerLimitClip + upperLimitClip)/2;

    cout << "should calc from " << lowerLimit << " to " << upperLimit << " ? " << endl;

    hist.save("/tmp/hist.arma");

    window.adapt(lowerLimit, upperLimit);

}

double System::getGlobalExtremum(const System::extrema type)
{
    using std::min_element;
    typedef std::function<bool(const double &, const double &)> compFuncType;


    double localExtrema, valueDifference;
    uint xd, yd, zd, particleIndex;

    xd = yd = zd = 0;

    //vector representing particles = 0, 1, ..., nParticles - 1
    vector<uint> particleIndices(m_nParticles);
    for (uint particleIndex = 0; particleIndex < m_nParticles; ++particleIndex)
    {
        particleIndices.at(particleIndex) = particleIndex;
    }

    compFuncType lessThan = [] (const double &v1, const double & v2) {return v1 < v2;};
    compFuncType greaterThan = [] (const double &v1, const double & v2) {return v1 > v2;};

    compFuncType arrangeSortCompare;
    compFuncType extremumCheckCompare;

    if (type == System::extrema::maximum)
    {
        extremumCheckCompare = greaterThan;
        arrangeSortCompare = lessThan;
    }
    else
    {
        extremumCheckCompare = lessThan;
        arrangeSortCompare = greaterThan;
    }

    localExtrema = 1.0;
    while (localExtrema != 0)
    {
        localExtrema = 0;

        //Arrange particles by values
        particleIndex = *(std::min_element(particleIndices.begin(),
                                           particleIndices.end(),
                                           [this, arrangeSortCompare] (const uint &p1, const uint &p2)
        {
            return arrangeSortCompare(getValue(p1), getValue(p2));
        }));


        //Search for the displacement which extremizes value gain
        for (uint x = 0; x < m_NX; ++x)
        {
            for (uint y = 0; y < m_NY; ++y)
            {
                for (uint z = 0; z < m_NZ; ++z)
                {
                    if (isOccupiedLoction(x, y, z))
                    {
                        continue;
                    }

                    valueDifference = getValueDifference(particleIndex, x, y, z);

                    if (extremumCheckCompare(valueDifference, localExtrema))
                    {
                        localExtrema = valueDifference;

                        xd = x;
                        yd = y;
                        zd = z;
                    }
                }
            }
        }

        //if no such displacement exist, we go on to the next particle.
        //if no such displacement exist for all particles, then we end.
        //if it does, we perform it, then reset the particle loop.
        if (localExtrema == 0)
        {
            break;
        }
        else
        {
            changePosition(particleIndex, xd, yd, zd);
        }
    }

    return getTotalValue();

}

void System::getRandomParticleAndDestination(uint &particleIndex, uint &xd, uint &yd, uint &zd) const
{
    particleIndex = m_URNG()*m_nParticles;
    uint destination = m_URNG()*m_freeVolume;

    findDestination(destination, xd, yd, zd);
}

void System::randomizeParticlePositions()
{
    uint xd, yd, zd, destination;

    xd = yd = zd = 0;

    uint nSweeps = 5;
    uint sweep = 0;

    while (sweep < nSweeps)
    {
        for (uint particleIndex = 0; particleIndex < m_nParticles; ++particleIndex)
        {

            destination = m_URNG()*m_freeVolume;

            findDestination(destination, xd, yd, zd);

            changePosition(particleIndex, xd, yd, zd);

        }

        sweep++;

    }

}

