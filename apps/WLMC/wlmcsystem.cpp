#include "wlmcsystem.h"
#include "wlmcwindow.h"

using namespace WLMC;

WLMCSystem::WLMCSystem(const uint nParticles,
                       const uint NX,
                       const uint NY,
                       const uint NZ,
                       const uint movesPerSampling,
                       const double flatnessCriterion,
                       const uint overlap,
                       const uint nbinsOverMinWindowSizeFlat,
                       const uint minWindowSizeRough,
                       const uint windowIncrementSize,
                       const double *f,
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
    m_nbinsOverMinWindowSizeFlat(nbinsOverMinWindowSizeFlat),
    m_minWindowSizeRough(minWindowSizeRough),
    m_windowIncrementSize(windowIncrementSize),
    m_f(f),
    m_URNG(URNG)
{

}

void WLMCSystem::sampleWindow(WLMCWindow *window)
{
    uint nMoves = 0;
    
    while (nMoves < m_movesPerSampling)
    {
        if (doSingleMove(window))
        {
            nMoves++;
        }
    }
}

bool WLMCSystem::doSingleMove(WLMCWindow *window)
{
    uint particleIndex, xd, yd, zd;

    getRandomParticleAndDestination(particleIndex, xd, yd, zd);

    double oldValue = getTotalValue();
    double newValue = oldValue + getValueDifference(particleIndex, xd, yd, zd);

    if (!window->isLegal(oldValue) || !window->isLegal(newValue))
    {
        changePosition(particleIndex, xd, yd, zd);

        return false;
    }

    uint oldBin = window->getBin(oldValue);
    uint newBin = window->getBin(newValue);

    double oldDOS = window->DOS(oldBin);
    double newDOS = window->DOS(newBin);

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

void WLMCSystem::findDestination(const uint destination, uint &xd, uint &yd, uint &zd)
{
    uint search = 0;

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

void WLMCSystem::locateGlobalExtremaValues(double &min, double &max, kMC::KMCSolver *solver)
{
    uint nSweeps = 1;
    uint sweep = 0;

    double localMax, localMin;

    min = std::numeric_limits<double>::max();
    max = std::numeric_limits<double>::min();

    vector<double> allExtrema;
    bool isIn;

    while (sweep < nSweeps)
    {
        randomizeParticlePositions();

        localMin = getGlobalExtremum(WLMCSystem::extrema::minimum);

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
            solver->dumpLAMMPS(allExtrema.size());
        }

        if (localMin < min)
        {
            min = localMin;
        }

        sweep++;
    }

    nSweeps = 100;
    sweep = 0;

    while (sweep < nSweeps)
    {
        randomizeParticlePositions();

        localMax = getGlobalExtremum(WLMCSystem::extrema::maximum);

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
            solver->dumpLAMMPS(allExtrema.size());
        }

        if (localMax > max)
        {
            max = localMax;
        }

        sweep++;
    }

    cout << "found " << allExtrema.size() << " extrema." << endl;

}

double WLMCSystem::getGlobalExtremum(const WLMCSystem::extrema type)
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

    if (type == WLMCSystem::extrema::maximum)
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

void WLMCSystem::getRandomParticleAndDestination(uint &particleIndex, uint &xd, uint &yd, uint &zd)
{
    particleIndex = m_URNG()*m_nParticles;
    uint destination = m_URNG()*m_freeVolume;

    findDestination(destination, xd, yd, zd);
}

void WLMCSystem::randomizeParticlePositions()
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

