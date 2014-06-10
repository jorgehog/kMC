#include "WLMCEvent.h"

using namespace kMC;



void WLMCEvent::initialize()
{

    m_minEnergies.set_size(solver()->volume());
    m_maxEnergies.set_size(solver()->volume());

    findAllExtrema();

    m_nCount = m_nStart;

    m_visitCounts.set_size(m_nbins);
    m_DOS.set_size(m_nbins);

    initializeNewCycle();


}

void WLMCEvent::execute()
{

    moveParticle();

    if (nTimesExecuted() % 1000 != 0 || nTimesExecuted() == 0)
    {
        return;
    }

    m_DOS = normalise(m_DOS);
    cout << "N = " << SoluteParticle::nParticles() << " f=" << m_f << endl;

    double flatness = estimateFlatness();

    setValue(flatness);

    output();

    if (flatness >= m_flatnessCriteria)
    {
        m_f = m_fIteratorFunction(m_f);
        resetCounts();

        initRandomConfiguration();

        if (m_f < m_fEnd)
        {
            prepNextOccupancyLevel();
        }
    }
}

void WLMCEvent::output() const
{
    stringstream file;
    file << solver()->filePath() << "/stateDensity" << SoluteParticle::nParticles() << ".arma";
    vec E = linspace<vec>(m_minEnergies(m_nCount), m_maxEnergies(m_nCount), m_nbins);

    uvec indices = find(m_visitCounts != m_unsetCount);
    uvec vc = m_visitCounts(indices);

    cout << find(m_visitCounts == m_unsetCount).t() << endl;

    vec dos = m_DOS(indices);
    vec e = E(indices);
    vec idx = conv_to<vec>::from(indices);
    vec vcd = conv_to<vec>::from(vc);

    if (join_rows(join_rows(join_rows(e, dos), vcd), idx).eval().save(file.str()))
    {
        cout << "storing " << file.str() << endl;
    }

}

double WLMCEvent::estimateFlatness() const
{

    uint min = arma::min(m_visitCounts);

    double mean = 0;
    uint c = 0;
    for (const uint &vc : m_visitCounts)
    {
        if (vc != m_unsetCount)
        {
            mean += vc;

            c++;
        }
    }

    if (c == 0)
    {
        cout << "no visited bins." << endl;
        return 0;
    }

    mean /= c;

    double flatness = min/mean;

    return flatness;
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
    if (SoluteParticle::nParticles() == NX()*NY()*NZ() - m_nSkipped + 2)
    {
        solver()->exit();
    }

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



