#include "kmcsolver.h"

#include "reactions/reaction.h"
#include "reactions/diffusion/diffusionreaction.h"

#include "boundary/boundary.h"

#include <sys/time.h>

#include <armadillo>

#include <iostream>
#include <iomanip>
#include <fstream>


using namespace arma;
using namespace std;
using namespace kMC;

KMCSolver::KMCSolver(const Setting & root)
{

    onConstruct();

    const Setting & SystemSettings = getSurfaceSetting(root, "System");
    const Setting & SolverSettings = getSurfaceSetting(root, "Solver");
    const Setting & diffusionSettings = getSetting(root, {"Reactions", "Diffusion"});


    Reaction::loadConfig(getSurfaceSetting(root, "Reactions"));

    DiffusionReaction::loadConfig(diffusionSettings);

    Site::loadConfig(SystemSettings);


    setNumberOfCycles(
                getSurfaceSetting<uint>(SolverSettings, "nCycles"));

    setCyclesPerOutput(
                getSurfaceSetting<uint>(SolverSettings, "cyclesPerOutput"));

    setRNGSeed(
                getSurfaceSetting<uint>(SolverSettings, "seedType"),
                getSurfaceSetting<int>(SolverSettings, "specificSeed"));

    setTargetSaturation(
                getSurfaceSetting<double>(SystemSettings, "SaturationLevel"));


    uvec3 boxSize;

    boxSize(0) = getSurfaceSetting(SystemSettings, "BoxSize")[0];
    boxSize(1) = getSurfaceSetting(SystemSettings, "BoxSize")[1];
    boxSize(2) = getSurfaceSetting(SystemSettings, "BoxSize")[2];

    setBoxSize(boxSize);

}

KMCSolver::KMCSolver()
{
    onConstruct();
}

KMCSolver::~KMCSolver()
{

    finalizeObject();

    clearSites();

    Site::clearAll();
    DiffusionReaction::clearAll();
    Boundary::clearAll();

    refCounter--;

}

void KMCSolver::checkRefCounter()
{
    if (refCounter > 1)
    {
        cerr << "ERROR: "<< refCounter << " solver objects alive when freeing.\n";
        cerr << "Static member variables of objects IN USE by living solver WILL BE FREED." << endl;

        exit();
    }
}

void KMCSolver::onConstruct()
{

    m_NX = UNSET_UINT;
    m_NY = UNSET_UINT;
    m_NZ = UNSET_UINT;

    outputCounter = 0;

    totalTime = 0;

    cycle = 1;

    m_kTot = 0;

    Boundary::setMainSolver(this);

    Reaction::setMainSolver(this);

    Site::setMainSolver(this);

    refCounter++;

}

void KMCSolver::finalizeObject()
{
    checkRefCounter();

    KMCDebugger_Finalize();
}

void KMCSolver::mainloop()
{

    Reaction * selectedReaction;
    uint choice;
    double R;

    dumpXYZ();

    KMCDebugger_Init();

    while(cycle <= m_nCycles)
    {

        getRateVariables();

        R = m_kTot*KMC_RNG_UNIFORM();

        choice = getReactionChoice(R);

        selectedReaction = m_allPossibleReactions.at(choice);
        KMCDebugger_SetActiveReaction(selectedReaction);

        selectedReaction->execute();
        KMCDebugger_PushTraces();


        if (cycle%m_cyclesPerOutput == 0)
        {
            dumpOutput();
            dumpXYZ();
        }


        totalTime += Reaction::linearRateScale()/m_kTot;
        cycle++;


        Site::updateBoundaries();

    }

}

void KMCSolver::reset()
{

    finalizeObject();

    //TMP
    m_kTot2 = 0;
    m_accuAllRates2.clear();
    m_allPossibleReactions2.clear();
    prevUpdatedReacs.clear();
    prevUpdatedReacsSet.clear();
    //

    totalTime = 0;

    cycle = 1;

    outputCounter = 0;

    m_kTot = 0;

    m_allPossibleReactions.clear();

    m_accuAllRates.clear();

    Site::clearAffectedSites();

    forEachSiteDo([] (Site * site)
    {
        site->reset();
    });

    KMCDebugger_Assert(accu(Site::totalActiveParticlesVector()), ==, 0);

    KMCDebugger_Assert(Site::totalDeactiveParticles(ParticleStates::solution), ==, m_NX*m_NY*m_NZ);

    KMCDebugger_AssertClose(Site::totalEnergy(), 0, 1E-5);

    KMCDebugger_Assert(Site::totalActiveSites(), ==, 0);

    Site::setZeroTotalEnergy();

    Site::finalizeBoundaries();

    setRNGSeed(Seed::specific, Seed::initialSeed);

    Site::initializeBoundaries();

    KMCDebugger_Init();

}



void KMCSolver::dumpXYZ()
{

    stringstream s;
    s << "kMC" << outputCounter++ << ".xyz";

    ofstream o;
    o.open("outfiles/" + s.str());



    stringstream surface;
    stringstream crystal;
    stringstream solution;

    uint nLines = 0;
    s.str(string());

    Site* currentSite;

    for (uint i = 0; i < m_NX; ++i)
    {
        for (uint j = 0; j < m_NY; ++j)
        {
            for (uint k = 0; k < m_NZ; ++k)
            {

                currentSite = sites[i][j][k];

                bool isSurface = currentSite->isSurface();

                if (currentSite->isActive() || isSurface)
                {
                    s << "\n" << ParticleStates::shortNames.at(sites[i][j][k]->particleState()) << " " << i << " " << j << " " << k << " " << sites[i][j][k]->nNeighbors();

                    if (isSurface)
                    {
                        surface << s.str();
                    }

                    else if (currentSite->isCrystal())
                    {
                        crystal << s.str();
                    }

                    else
                    {
                        solution << s.str();
                    }

                    s.str(string());
                    nLines++;
                }
            }
        }
    }

    o << nLines << "\n - " << surface.str() << crystal.str() << solution.str();
    o.close();

}

void KMCSolver::dumpOutput()
{

    cout << setw(5) << right << setprecision(1) << fixed
         << (double)cycle/m_nCycles*100 << "%   "
         << outputCounter
         << endl;
    cout << setprecision(6);
}


void KMCSolver::forEachSiteDo(function<void (Site *)> applyFunction) const
{
    for (uint x = 0; x < m_NX; ++x)
    {
        for (uint y = 0; y < m_NY; ++y)
        {
            for (uint z = 0; z < m_NZ; ++z)
            {
                applyFunction(sites[x][y][z]);
            }
        }
    }
}

void KMCSolver::forEachSiteDo_sendIndices(function<void (Site *, uint, uint, uint)> applyFunction) const
{
    for (uint x = 0; x < m_NX; ++x)
    {
        for (uint y = 0; y < m_NY; ++y)
        {
            for (uint z = 0; z < m_NZ; ++z)
            {
                applyFunction(sites[x][y][z], x, y, z);
            }
        }
    }
}

void KMCSolver::forEachActiveSiteDo(function<void (Site *)> applyFunction) const
{
    for (uint x = 0; x < m_NX; ++x)
    {
        for (uint y = 0; y < m_NY; ++y)
        {
            for (uint z = 0; z < m_NZ; ++z)
            {
                if (sites[x][y][z]->isActive())
                {
                    applyFunction(sites[x][y][z]);
                }
            }
        }
    }
}

void KMCSolver::forEachActiveSiteDo_sendIndices(function<void (Site *, uint, uint, uint)> applyFunction) const
{
    for (uint x = 0; x < m_NX; ++x)
    {
        for (uint y = 0; y < m_NY; ++y)
        {
            for (uint z = 0; z < m_NZ; ++z)
            {
                if (sites[x][y][z]->isActive())
                {
                    applyFunction(sites[x][y][z], x, y, z);
                }
            }
        }
    }
}


void KMCSolver::initializeSites()
{

    sites = new Site***[m_NX];

    for (uint x = 0; x < m_NX; ++x)
    {
        sites[x] = new Site**[m_NY];

        for (uint y = 0; y < m_NY; ++y)
        {
            sites[x][y] = new Site*[m_NZ];

            for (uint z = 0; z < m_NZ; ++z)
            {
                sites[x][y][z] = new Site(x, y, z);
            }
        }
    }


    initializeSiteNeighborhoods();

}


void KMCSolver::clearSites()
{

    KMCDebugger_SetEnabledTo(false);

    for (uint i = 0; i < m_NX; ++i)
    {
        for (uint j = 0; j < m_NY; ++j)
        {
            for (uint k = 0; k < m_NZ; ++k)
            {
                delete sites[i][j][k];
            }

            delete [] sites[i][j];
        }

        delete [] sites[i];
    }

    delete [] sites;


    KMCDebugger_Assert(accu(Site::totalActiveParticlesVector()), ==, 0);

    KMCDebugger_Assert(accu(Site::totalDeactiveParticlesVector()), ==, 0);

    KMCDebugger_AssertClose(Site::totalEnergy(), 0, 1E-5);

    KMCDebugger_Assert(Site::totalActiveSites(), ==, 0);

    Site::clearAffectedSites();
    Site::setZeroTotalEnergy();

    Reaction::clearAll();

    m_allPossibleReactions.clear();


    KMCDebugger_ResetEnabled();

}

void KMCSolver::setBoxSize_KeepSites(const uvec3 &boxSizes)
{

    // have to clear sites and all that jazz because of vector setup...

}



void KMCSolver::initializeCrystal(const double relativeSeedSize)
{

    if (relativeSeedSize > 1.0)
    {
        cerr << "The seed size cannot exceed the box size." << endl;
        KMCSolver::exit();
    }

    else if (relativeSeedSize < 0)
    {
        cerr << "The seed size cannot be negative." << endl;
        KMCSolver::exit();
    }

    bool noSeed = false;
    KMCDebugger_SetEnabledTo(false);

    if (!noSeed)
    {
        sites[m_NX/2][m_NY/2][m_NZ/2]->spawnAsFixedCrystal();
        KMCDebugger_PushTraces();
    }

    uint crystalSizeX = round(m_NX*relativeSeedSize);
    uint crystalSizeY = round(m_NY*relativeSeedSize);
    uint crystalSizeZ = round(m_NZ*relativeSeedSize);

    uint crystalStartX = m_NX/2 - crystalSizeX/2;
    uint crystalStartY = m_NY/2 - crystalSizeY/2;
    uint crystalStartZ = m_NZ/2 - crystalSizeZ/2;

    uint crystalEndX = crystalStartX + crystalSizeX;
    uint crystalEndY = crystalStartY + crystalSizeY;
    uint crystalEndZ = crystalStartZ + crystalSizeZ;

    uint solutionEndX = crystalStartX - Site::nNeighborsLimit();
    uint solutionEndY = crystalStartY - Site::nNeighborsLimit();
    uint solutionEndZ = crystalStartZ - Site::nNeighborsLimit();

    uint solutionStartX = crystalEndX + Site::nNeighborsLimit();
    uint solutionStartY = crystalEndY + Site::nNeighborsLimit();
    uint solutionStartZ = crystalEndZ + Site::nNeighborsLimit();



    for (uint i = 0; i < m_NX; ++i)
    {
        for (uint j = 0; j < m_NY; ++j)
        {
            for (uint k = 0; k < m_NZ; ++k)
            {

                if (!noSeed)
                {

                    if (i >= crystalStartX && i < crystalEndX)
                    {
                        if (j >= crystalStartY && j < crystalEndY)
                        {
                            if (k >= crystalStartZ && k < crystalEndZ)
                            {
                                if (!((i == m_NX/2 && j == m_NY/2 && k == m_NZ/2)))
                                {
                                    sites[i][j][k]->activate();
                                    KMCDebugger_PushTraces();
                                }

                                continue;

                            }
                        }
                    }

                }

                if (noSeed || ((i < solutionEndX) || (i >= solutionStartX) || (j < solutionEndY) || (j >= solutionStartY) || (k < solutionEndZ) || (k >= solutionStartZ)))
                {
                    if (KMC_RNG_UNIFORM() < m_targetSaturation)
                    {
                        if(sites[i][j][k]->isLegalToSpawn())
                        {
                            sites[i][j][k]->activate();
                            KMCDebugger_PushTraces();
                        }
                    }
                }
            }
        }
    }

    KMCDebugger_ResetEnabled();

}

void KMCSolver::initializeSolutionBath()
{

    forEachSiteDo([this] (Site * site)
    {
        if (site->isLegalToSpawn())
        {
            if (KMC_RNG_UNIFORM() < targetSaturation())
            {
                site->activate();
            }
        }
    });
}




void KMCSolver::getRateVariables()
{

    m_kTot = 0;

    m_accuAllRates.clear();

    m_allPossibleReactions.clear();


    Site::updateAffectedSites();

    double minRate = std::numeric_limits<double>::max();
    forEachSiteDo([this, &minRate] (Site * site)
    {
        site->forEachActiveReactionDo([this, &minRate] (Reaction * reaction)
        {

            KMCDebugger_Assert(reaction->rate(), !=, Reaction::UNSET_RATE, "Reaction rate should not be unset at this point.", reaction->getFinalizingDebugMessage());

            m_kTot += reaction->rate();

            if (reaction->rate() < minRate)
            {
                minRate = reaction->rate();
            }

            m_accuAllRates.push_back(m_kTot);

            m_allPossibleReactions.push_back(reaction);

//            KMCDebugger_Assert(m_accuAllRates.at(m_accuAllRates.size()-1), ==, m_accuAllRates2.at(m_accuAllRates.size()-1));
//            KMCDebugger_Assert(m_allPossibleReactions.at(m_allPossibleReactions.size()-1), ==, reaction);

        });
    });

    KMCDebugger_AssertClose(m_kTot, m_kTot2, minRate/2);

}


uint KMCSolver::getReactionChoice(double R)
{

    KMCDebugger_Assert(m_accuAllRates.size(), !=, 0, "No active reactions.");

    uint imax = m_accuAllRates.size() - 1;
    uint MAX = imax;
    uint imin = 0;
    uint imid = 1;

    // continue searching while imax != imin + 1
    while (imid != imin)
    {

        // calculate the midpoint for roughly equal partition
        imid = imin + (imax - imin)/2;

        //Is the upper limit above mid?
        if (R > m_accuAllRates.at(imid))
        {

            if (imid == MAX)
            {
                return MAX;
            }

            //Are we just infront of the limit?
            else if (R < m_accuAllRates.at(imid + 1))
            {
                //yes we were! Returning current mid + 1.
                //If item i in accuAllrates > R, then reaction i is selected.
                //This because there is no zero at the start of accuAllrates.

                return imid + 1;
            }

            //No we're not there yet, so we search above us.
            else
            {
                imin = imid + 1;
            }
        }

        //No it's not. Starting new search below mid!
        else
        {

            if (imid == 0)
            {
                return 0;
            }

            imax = imid;
        }


    }

    return imid + 1;

}

void KMCSolver::setBoxSize(const uvec3 boxSize, bool check, bool keepSystem)
{

    if (keepSystem)
    {
        setBoxSize_KeepSites(boxSize);
        return;
    }

    if (m_NX != UNSET_UINT && m_NY != UNSET_UINT && m_NZ != UNSET_UINT)
    {
        clearSites();
    }

    m_NX = boxSize(0);
    m_NY = boxSize(1);
    m_NZ = boxSize(2);


    m_N = boxSize;

    if (Site::nNeighborsLimit() != UNSET_UINT && check)
    {
        if (Site::nNeighborsLimit() >= min(m_N)/2)
        {
            cerr << "Neighbor reach must be lower than half the minimum box dimension to avoid sites directly affecting themselves." << endl;
            KMCSolver::exit();
        }
    }


    initializeSites();

    Site::initializeBoundaries();

    initializeDiffusionReactions();

}

void KMCSolver::setRNGSeed(uint seedState, int defaultSeed)
{

    seed_type prevSeed = Seed::initialSeed;

    seed_type seed = -1;

    switch (seedState)
    {
    case Seed::specific:
        seed = static_cast<seed_type>(defaultSeed);
        break;
    case Seed::fromTime:
        seed = static_cast<seed_type>(time(NULL));
        break;
    default:
        throw std::runtime_error("Seed not specified.");
        break;
    }

    KMC_INIT_RNG(seed);


    if (prevSeed != Seed::initialSeed)
    {
        cout << " -- new seed set: " << seed << endl;
    }

}


uint KMCSolver::refCounter = 0;
