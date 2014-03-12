#include "kmcsolver.h"

#include "reactions/reaction.h"
#include "reactions/diffusion/diffusionreaction.h"

#include "boundary/boundary.h"

#include "debugger/debugger.h"

#include <sys/time.h>

#include <armadillo>

#include <iostream>
#include <iomanip>
#include <fstream>


using namespace arma;
using namespace std;
using namespace kMC;

KMCSolver::KMCSolver(const Setting & root) :
    m_NX(UNSET_UINT),
    m_NY(UNSET_UINT),
    m_NZ(UNSET_UINT),
    totalTime(0),
    cycle(0),
    outputCounter(0)
{

    const Setting & SystemSettings = getSurfaceSetting(root, "System");
    const Setting & SolverSettings = getSurfaceSetting(root, "Solver");
    const Setting & InitializationSettings = getSurfaceSetting(root, "Initialization");
    const Setting & diffusionSettings = getSetting(root, {"Reactions", "Diffusion"});


    Boundary::setMainSolver(this);

    Reaction::setMainSolver(this);

    Site::setMainSolver(this);



    Reaction::loadConfig(getSurfaceSetting(root, "Reactions"));

    DiffusionReaction::loadConfig(diffusionSettings);

    Site::loadConfig(SystemSettings);


    setNumberOfCycles(
                getSurfaceSetting<uint>(SolverSettings, "nCycles"));

    setCyclesPerOutput(
                getSurfaceSetting<uint>(SolverSettings, "cyclesPerOutput"));

    setSaturation(
                getSurfaceSetting<double>(InitializationSettings, "SaturationLevel"));

    setRelativeSeedSize(
                getSurfaceSetting<double>(InitializationSettings, "RelativeSeedSize"));

    setRNGSeed(
                getSurfaceSetting<uint>(SolverSettings, "seedType"),
                getSurfaceSetting<int>(SolverSettings, "specificSeed"));

    Site::setInitialBoundaries(
                getSurfaceSetting(SystemSettings, "Boundaries"));

    DiffusionReaction::setSeparation(
                getSurfaceSetting<uint>(diffusionSettings, "separation"));



    uvec3 boxSize;

    boxSize(0) = getSurfaceSetting(SystemSettings, "BoxSize")[0];
    boxSize(1) = getSurfaceSetting(SystemSettings, "BoxSize")[1];
    boxSize(2) = getSurfaceSetting(SystemSettings, "BoxSize")[2];

    setBoxSize(boxSize);



    ptrCount++;

}

KMCSolver::~KMCSolver()
{

    if (ptrCount > 1)
    {
        cout << "WARNING: "<< ptrCount << " solver objects alive when freeing.\n";
        cout << "Static member variables of objects IN USE by living solver WILL BE FREED." << endl;
    }

    clearSites();

    m_allReactions.clear();

    Site::resetAll();
    Reaction::resetAll();
    DiffusionReaction::resetAll();
    Boundary::resetAll();

    KMCDebugger_Finalize();

    ptrCount--;

}



void KMCSolver::run()
{

    KMCDebugger_Init();

    Reaction * selectedReaction;
    uint choice;
    double R;

    initializeCrystal();

    while(cycle < m_nCycles)
    {

        getRateVariables();

        double r = KMC_RNG_UNIFORM();

        R = m_kTot*r;

        choice = getReactionChoice(R);

        selectedReaction = m_allReactions.at(choice);
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

    KMCDebugger_Init();

    clearSites();

    initializeSites();

    Site::initializeBoundaries();

    initializeDiffusionReactions();

    setRNGSeed(Seed::specific, Seed::initialSeed);


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


void KMCSolver::initializeDiffusionReactions()
{

    //Loop over all sites
    for (uint x = 0; x < m_NX; ++x)
    {
        for (uint y = 0; y < m_NY; ++y)
        {
            for (uint z = 0; z < m_NZ; ++z)
            {

                sites[x][y][z]->initializeDiffusionReactions();

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

void KMCSolver::clearSiteNeighborhoods()
{
    for (uint x = 0; x < m_NX; ++x)
    {
        for (uint y = 0; y < m_NY; ++y)
        {
            for (uint z = 0; z < m_NZ; ++z)
            {
                sites[x][y][z]->clearNeighborhood();
            }
        }
    }
}

void KMCSolver::initializeSiteNeighborhoods()
{
    for (uint x = 0; x < m_NX; ++x)
    {
        for (uint y = 0; y < m_NY; ++y)
        {
            for (uint z = 0; z < m_NZ; ++z)
            {
                sites[x][y][z]->introduceNeighborhood();
            }
        }
    }
}

void KMCSolver::clearSites()
{

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

    KMCDebugger_AssertClose(Site::totalEnergy(), 0, 1E-5);

    KMCDebugger_Assert(Site::totalActiveSites(), ==, 0);

    Site::resetAffectedSites();
    Site::setZeroTotalEnergy();

}


void KMCSolver::initializeCrystal()
{

    bool enabled = KMCDebugger_IsEnabled;
    KMCDebugger_SetEnabledTo(false);

    sites[m_NX/2][m_NY/2][m_NZ/2]->spawnAsFixedCrystal();
    KMCDebugger_PushTraces();

    uint crystalSizeX = round(m_NX*m_relativeSeedSize);
    uint crystalSizeY = round(m_NY*m_relativeSeedSize);
    uint crystalSizeZ = round(m_NZ*m_relativeSeedSize);

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

                if ((i < solutionEndX) || (i >= solutionStartX) || (j < solutionEndY) || (j >= solutionStartY) || (k < solutionEndZ) || (k >= solutionStartZ))
                {
                    if (KMC_RNG_UNIFORM() < m_saturation)
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


    dumpXYZ();


    KMCDebugger_SetEnabledTo(enabled);

}




void KMCSolver::getRateVariables()
{

    m_kTot = 0;
    m_accuAllRates.clear();
    m_allReactions.clear();

    Site::updateAffectedSites();

    for (uint x = 0; x < m_NX; ++x)
    {
        for (uint y = 0; y < m_NY; ++y)
        {
            for (uint z = 0; z < m_NZ; ++z)
            {
                for (Reaction* reaction : sites[x][y][z]->activeReactions())
                {
                    assert(reaction->rate() != Reaction::UNSET_RATE);
                    m_kTot += reaction->rate();
                    m_accuAllRates.push_back(m_kTot);
                    m_allReactions.push_back(reaction);
                }
            }
        }
    }

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

void KMCSolver::setBoxSize(const uvec3 boxSize)
{

    if (m_NX != UNSET_UINT && m_NY != UNSET_UINT && m_NZ != UNSET_UINT)
    {
        clearSites();
    }

    m_NX = boxSize(0);
    m_NY = boxSize(1);
    m_NZ = boxSize(2);


    m_N = boxSize;

    if (Site::nNeighborsLimit() != UNSET_UINT)
    {
        if (Site::nNeighborsLimit() >= min(m_N)/2)
        {
            cerr << "Neighbor reach must be lower than half the minimum box dimension to avoid sites directly affecting themselves." << endl;
            exit(1);
        }
    }

    initializeSites();

    Site::initializeBoundaries();

    initializeDiffusionReactions();

}

void KMCSolver::setRNGSeed(uint seedState, int defaultSeed = 0)
{

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
        cerr << "SEED NOT SPECIFIED." << endl;
        throw Seed::seedNotSetException;
        break;
    }

    //    cout << "initializing seed : " << seed << endl;
    KMC_INIT_RNG(seed);

}


uint KMCSolver::ptrCount = 0;
