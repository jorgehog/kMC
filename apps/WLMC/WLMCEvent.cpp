#include "WLMCEvent.h"

#include "wlmcwindow.h"
#include "kmcwlmcsystem.h"


using namespace kMC;
using namespace WLMC;


void WLMCEvent::initialize()
{

    m_minBin = 0;
    m_maxBin = m_nbins - 1;

    m_minEnergies.set_size(solver()->volume());
    m_maxEnergies.set_size(solver()->volume());

//    findAllExtrema();

    m_nCount = m_nStart;

    m_visitCounts.set_size(m_nbins);
    m_DOS.set_size(m_nbins);

//    initializeNewCycle();

}

void WLMCEvent::execute()
{
    initRandomConfiguration();

    double f = m_fBegin;

    KMCWLMCSystem system(solver(),
                         m_movesPerWindowCheck,
                         m_flatnessCriteria,
                         m_windowOverlap,
                         m_nbinsOverMinWindowSizeFlat,
                         m_minWindowSizeRough,
                         m_windowIncrementSize,
                         &f);

    double max, min;
    system.locateGlobalExtremaValues(min, max);

//    cout << min << " " << m_minEnergies(m_nCount) << endl;
//    cout << max << " " << m_maxEnergies(m_nCount) << endl;

//    min = 377.383;
//    max = 547.847;


    WLMCWindow mainWindow(&system, m_nbins, min, max);

    uint lowerClip, upperClip;
    system.clipWindow(lowerClip, upperClip, mainWindow);

    WLMCWindow cheatWindow(&system,
                           mainWindow.DOS(),
                           mainWindow.energies(),
                           lowerClip,
                           upperClip,
                           WLMCWindow::OVERLAPTYPES::NONE);

    system.setupPresetWindowConfigurations(cheatWindow.minValue(), cheatWindow.maxValue(), 10);

    while (f >= m_fEnd)
    {
        cheatWindow.calculateWindow();

        f = m_fIteratorFunction(f);

        cheatWindow.reset();
    }

    return;

    while (m_f > m_fEnd)
    {
        calculateWindow(m_minBin, m_maxBin);
        prepNextOccupancyLevel();
    }

}

void WLMCEvent::calculateWindow(const uint startBin, const uint endBin)
{
    m_windowMinBin = startBin;
    m_windowMaxBin = endBin;

    bool flat = false;

    uint counter = 0;

    uint flatSubwindowStart;
    uint flatSubwindowEnd;

    while (!flat)
    {
        counter++;

        moveParticle();
        setValue(estimateFlatness(m_visitCounts));

        if (counter % m_movesPerWindowCheck != 0)
        {
            continue;
        }


        m_DOS = normalise(m_DOS);

        findFlatWindow(m_visitCounts(span(m_windowMinBin, m_windowMaxBin)), flatSubwindowStart, flatSubwindowEnd, pointClosestToMean());

        flat = flatSubwindowEnd - flatSubwindowStart >= m_nbins/m_nbinsOverMinWindowSizeFlat && false;

        ofstream f;
        f.open(solver()->filePath() + "/flatness.txt");
        f << flatSubwindowStart << "\n" << flatSubwindowEnd << endl;
        f.close();

        output();
    }



}

uint WLMCEvent::pointClosestToMean() const
{
    double mean = 0;
    uint count = 0;

    for (uint i = m_windowMinBin; i <= m_windowMaxBin; ++i)
    {
        if (m_visitCounts(i) == m_unsetCount)
        {
            continue;
        }

        mean += m_visitCounts(i);
        count++;

    }

    uint point = 0;
    double min = numeric_limits<double>::max();

    for (uint i = m_windowMinBin; i <= m_windowMaxBin; ++i)
    {
        if (m_visitCounts(i) == m_unsetCount)
        {
            continue;
        }

        double div = (m_visitCounts(i) - mean)*(m_visitCounts(i) - mean);

        if (div < min)
        {
            min = div;
            point = i;
        }

    }

    return point;
}

void WLMCEvent::output() const
{
    stringstream file;
    file << solver()->filePath() << "/stateDensity" << SoluteParticle::nParticles() << ".arma";
    vec E = linspace<vec>(m_minEnergies(m_nCount), m_maxEnergies(m_nCount), m_nbins);

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

double WLMCEvent::estimateFlatness(const uvec &visitCounts) const
{

    uint min = arma::min(visitCounts);

    double mean = 0;
    uint c = 0;
    for (const uint &vc : visitCounts)
    {
        if (vc != m_unsetCount)
        {
            mean += vc;

            c++;
        }
    }

    if (c == 0)
    {
        return 0;
    }

    mean /= c;

    double flatness = min/mean;

    return flatness;
}

bool WLMCEvent::isFlat(const uvec &visitCounts) const
{
    return estimateFlatness(visitCounts) > m_flatnessCriteria;
}

void WLMCEvent::findFlatWindow(const uvec &visitCounts, uint &lowerLimit, uint &upperLimit, const uint origin) const
{
    uint top = visitCounts.n_elem - 1;
    uint bottom = 0;

    if (origin > top - (m_nbins/m_nbinsOverMinWindowSizeFlat)/2)
    {
        upperLimit = top;
        lowerLimit = origin - m_nbinsOverMinWindowSizeFlat;
    }

    else if (origin < bottom + (m_nbins/m_nbinsOverMinWindowSizeFlat)/2)
    {
        lowerLimit = bottom;
        upperLimit = bottom + (m_nbins/m_nbinsOverMinWindowSizeFlat);
    }

    else
    {
        upperLimit = origin + (m_nbins/m_nbinsOverMinWindowSizeFlat)/2;
        lowerLimit = origin - (m_nbins/m_nbinsOverMinWindowSizeFlat)/2;
    }

    cout << origin << " " << lowerLimit << " " << upperLimit << endl;

    cout << "- " << endl;

    bool upperSet = false;
    bool lowerSet = false;

    //right propagation
    while (!(upperSet && lowerSet))
    {

        KMCDebugger_Assert(lowerLimit, <, upperLimit);

        if (!upperSet)
        {

            if (upperLimit == top)
            {
                upperSet = true;
            }

            if (!isFlat(visitCounts(span(lowerLimit, topIncrement(upperLimit, top)))))
            {
                upperSet = true;
            }
            else
            {
                upperLimit = topIncrement(upperLimit, top);
            }

        }

        if (!lowerSet)
        {

            if (lowerLimit == bottom)
            {
                lowerSet = true;
            }

            if (!isFlat(visitCounts(span(bottomIncrement(lowerLimit, bottom), upperLimit))))
            {
                lowerSet = true;
            }
            else
            {
                lowerLimit = bottomIncrement(lowerLimit, bottom);
            }

        }

        cout << lowerLimit << " " << upperLimit << endl;


    }

}


void WLMCEvent::initializeNewCycle()
{
    resetCounts();

    m_DOS.ones();
    m_f = m_fBegin;

    m_energySpan = m_maxEnergies(m_nCount) - m_minEnergies(m_nCount);

    initRandomConfiguration();

    KMCDebugger_Assert(m_energySpan, !=, 0);
}

void WLMCEvent::initRandomConfiguration()
{
    solver()->clearParticles();

    for (uint c = 0; c < m_nCount; ++c)
    {
        solver()->insertRandomParticle(0, false, false);
    }
}

void WLMCEvent::moveParticle()
{
    uint which = KMC_RNG_UNIFORM()*SoluteParticle::nParticles();
    uint where = KMC_RNG_UNIFORM()*(solver()->volume() - SoluteParticle::nParticles());

    SoluteParticle* particle = solver()->particle(which);

    uint xd = 0;
    uint yd = 0;
    uint zd = 0;

    findAvailableSite(where, xd, yd, zd);

    performTrialMove(particle, xd, yd, zd);

}

void WLMCEvent::findAvailableSite(const uint where, uint &xd, uint &yd, uint &zd) const
{
    uint search = 0;

    for (uint x = 0; x < NX(); ++x)
    {
        for (uint y = 0; y < NY(); ++y)
        {
            for (uint z = 0; z < NZ(); ++z)
            {
                if (!solver()->getSite(x, y, z)->isActive())
                {
                    if (search == where)
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

void WLMCEvent::performTrialMove(SoluteParticle *particle, const uint xd, const uint yd, const uint zd)
{

    double eNew = SoluteParticle::totalEnergy() + getEnergyDifference(particle, xd, yd, zd);

    uint prevBin = getBin(SoluteParticle::totalEnergy());
    uint newBin = getBin(eNew);

    double prevDOS = m_DOS(prevBin);
    double newDOS = m_DOS(newBin);

    if (!isLegalBin(newBin) || !isLegalBin(prevBin))
    {
        particle->changePosition(xd, yd, zd);
        return;
    }

    bool accepted = true;

    if (prevDOS < newDOS)
    {
        accepted = (KMC_RNG_UNIFORM() < prevDOS/newDOS);
    }

    if (accepted)
    {
        particle->changePosition(xd, yd, zd);

        KMCDebugger_AssertClose(SoluteParticle::totalEnergy(), eNew, 1E-3);

        registerVisit(newBin);
    }

    else
    {
        registerVisit(prevBin);
    }

}

void WLMCEvent::registerVisit(const uint bin)
{

    if (m_visitCounts(bin) == m_unsetCount)
    {
        m_visitCounts(bin) = 1;
    }
    else
    {
        m_visitCounts(bin)++;
    }

    m_DOS(bin) *= m_f;


    KMCDebugger_Assert(m_visitCounts(bin), !=, m_unsetCount);
}

uint WLMCEvent::getBin(const double energy)
{
    uint bin = m_nbins*(energy - m_minEnergies(m_nCount))/m_energySpan;

    if (bin == m_nbins)
    {
        KMCDebugger_AssertClose(energy, m_maxEnergies(m_nCount), 1E-3);

        return m_nbins - 1;
    }
    else if (bin > m_nbins)
    {
        solver()->exit("illegal bin");
    }

    return bin;
}

void WLMCEvent::prepNextOccupancyLevel()
{

    m_f = m_fIteratorFunction(m_f);

    findNewEnergyExtrema();

    bool spawned = false;
    solver()->forEachSiteDo([&spawned, this] (uint x, uint y, uint z, Site *site)
    {
        if (!site->isActive() && !spawned)
        {
            spawned = true;
            solver()->forceSpawnParticle(x, y, z);
        }
    });

    m_nCount++;

    initializeNewCycle();

}

void WLMCEvent::parseForExtrema(const int sign)
{

    double localExtrema;
    uint xd, yd, zd;

    bool done;

    do {

        done = true;

        if (sign == -1)
        {
            solver()->sortParticles([] (SoluteParticle *p1, SoluteParticle *p2) {return p1->energy() > p2->energy();});
        }
        else
        {
            solver()->sortParticles([] (SoluteParticle *p1, SoluteParticle *p2) {return p1->energy() < p2->energy();});
        }

        for (SoluteParticle *particle : solver()->particles())
        {
            localExtrema = 0;

            solver()->forEachSiteDo([&] (uint x, uint y, uint z, Site* site)
            {
                if (site->isActive())
                {
                    return;
                }

                double dE = getEnergyDifference(particle, x, y, z);

                if (sign == 1 ? dE > localExtrema : dE < localExtrema)
                {
                    localExtrema = dE;

                    xd = x;
                    yd = y;
                    zd = z;

                }

            });

            if (localExtrema == 0)
            {
                continue;
            }

            done = false;

            particle->changePosition(xd, yd, zd);

        }

    } while (!done);


    solver()->dumpLAMMPS(2*(SoluteParticle::nParticles() - 2) + (1 + sign)/2);

}

double WLMCEvent::getEnergyDifference(const SoluteParticle *particle, const uint xd, const uint yd, const uint zd) const
{

    double eNew = 0;

    Site::forEachNeighborDo_sendIndices(xd, yd, zd, [&particle, &eNew] (Site *neighbor, uint i, uint j, uint k)
    {
        if (neighbor->isActive())
        {
            if (neighbor->associatedParticle() != particle)
            {
                eNew += DiffusionReaction::potential(i, j, k);
            }
        }
    });

    return 2*(eNew - particle->energy());
}

void WLMCEvent::findAllExtrema()
{
    uint n;

    solver()->forceSpawnParticle(0, 0, 0);
    solver()->forceSpawnParticle(NX() - 1, NY() - 1, NZ() - 1);

    n = 2;
    do
    {
        parseForExtrema(-1);
        m_minEnergies(n) = SoluteParticle::totalEnergy();

        n++;

        solver()->insertRandomParticle(0, false, false);

    } while (n < m_minEnergies.n_elem - m_nSkipped);

    solver()->clearParticles();


    solver()->forceSpawnParticle(NX()/2, NY()/2, NZ()/2);
    solver()->forceSpawnParticle(NX()/2 + 1, NY()/2, NZ()/2);

    n = 2;
    do
    {
        parseForExtrema(+1);
        m_maxEnergies(n) = SoluteParticle::totalEnergy();

        n++;

        solver()->insertRandomParticle(0, false, false);

    } while (n < m_maxEnergies.n_elem);

    solver()->clearParticles();
}

void WLMCEvent::findNewEnergyExtrema()
{
    double treshold = 1E-10;

    m_minBin = 0;
    while (m_DOS(m_minBin) < treshold)
    {
        m_minBin++;

        KMCDebugger_Assert(m_minBin, <, m_nbins);
    }

    m_maxBin = m_nbins - 1;
    while(m_DOS(m_maxBin) < treshold)
    {
        m_maxBin--;

        KMCDebugger_Assert(m_maxBin, >, 0);
    }

}

bool WLMCEvent::isLegalBin(const uint bin)
{
    return bin >= m_windowMinBin && bin <= m_windowMaxBin;
}



