#include "kmcsolver.h"

#include "RNG/kMCRNG.h"

#include "site.h"

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
    totalTime(0),
    cycle(0),
    outputCounter(0)
{

    const Setting & SystemSettings = getSurfaceSetting(root, "System");
    const Setting & SolverSettings = getSurfaceSetting(root, "Solver");
    const Setting & InitializationSettings = getSurfaceSetting(root, "Initialization");

    const Setting & BoxSize = getSurfaceSetting(SystemSettings, "BoxSize");


    NX = BoxSize[0];
    NY = BoxSize[1];
    NZ = BoxSize[2];

    Boundary::setMainSolver(this);

    Site::setMainSolver(this);
    Site::loadConfig(SystemSettings);

    Reaction::setMainSolver(this);
    Reaction::loadConfig(getSurfaceSetting(root, "Reactions"));

    DiffusionReaction::loadConfig(getSetting(root, {"Reactions", "Diffusion"}));

    m_nCycles = getSurfaceSetting<uint>(SolverSettings, "nCycles");

    cyclesPerOutput = getSurfaceSetting<uint>(SolverSettings, "cyclesPerOutput");


    Seed::SeedState seedState = static_cast<Seed::SeedState>(getSurfaceSetting<uint>(SolverSettings, "seedType"));

    seed_type seed = -1;
    switch (seedState)
    {
    case Seed::specific:
        seed = getSurfaceSetting<seed_type>(SolverSettings, "specificSeed");
        break;
    case Seed::fromTime:
        seed = static_cast<seed_type>(time(NULL));
        break;
    default:
        cout << "SEED NOT SPECIFIED." << endl;
        throw Seed::seedNotSetException;
        break;
    }

    cout << "initializing seed : " << seed << endl;
    KMC_INIT_RNG(seed);


    saturation = getSurfaceSetting<double>(InitializationSettings, "SaturationLevel");
    RelativeSeedSize = getSurfaceSetting<double>(InitializationSettings, "RelativeSeedSize");

    KMCDebugger_Assert(RelativeSeedSize, <, 1.0, "The seed size cannot exceed the box size.");


    initializeSites();

    initializeDiffusionReactions();

    KMCDebugger_Init();
    Site::initializeBoundaries();

    ptrCount++;

}

KMCSolver::~KMCSolver()
{

    if (ptrCount > 1)
    {
        cout << "WARNING: Several solver objects alive when freeing.";
        cout << "Static member variables of objects IN USE by living solver WILL BE FREED." << endl;
    }

    for (uint i = 0; i < NX; ++i)
    {
        for (uint j = 0; j < NY; ++j)
        {
            for (uint k = 0; k < NZ; ++k)
            {
                delete sites[i][j][k];
            }

            delete [] sites[i][j];
        }

        delete [] sites[i];
    }

    delete [] sites;

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


        if (cycle%cyclesPerOutput == 0)
        {
            dumpOutput();
            dumpXYZ();
        }


        totalTime += Reaction::linearRateScale()/m_kTot;
        cycle++;


        Site::updateBoundaries();

    }

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

    for (uint i = 0; i < NX; ++i)
    {
        for (uint j = 0; j < NY; ++j)
        {
            for (uint k = 0; k < NZ; ++k)
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

    Site* currentSite;
    Site* destination;

    //Loop over all sites
    for (uint x = 0; x < NX; ++x)
    {
        for (uint y = 0; y < NY; ++y)
        {
            for (uint z = 0; z < NZ; ++z)
            {

                currentSite = sites[x][y][z];

                assert(currentSite->siteReactions().size() == 0 && "Sitereactions are already set");

                //For each site, loop over all neightbours
                for (uint i = 0; i < 3; ++i)
                {
                    for (uint j = 0; j < 3; ++j)
                    {
                        for (uint k = 0; k < 3; ++k)
                        {

                            destination = currentSite->neighborHood(Site::nNeighborsLimit() - 1 + i,
                                                                    Site::nNeighborsLimit() - 1 + j,
                                                                    Site::nNeighborsLimit() - 1 + k);

                            //This means that the destination is blocked by boundaries
                            if (destination != NULL)
                            {
                                //This menas we are not at the current site.
                                if(destination != currentSite)
                                {
                                    DiffusionReaction* diffusionReaction = new DiffusionReaction(currentSite, destination);
                                    currentSite->addReaction(diffusionReaction);
                                }

                                else
                                {
                                    assert((i == 1) && (j == 1) && (k == 1));
                                }
                            }

                        }
                    }
                }
            }
        }
    }

}

void KMCSolver::initializeSites()
{

    sites = new Site***[NX];

    for (uint x = 0; x < NX; ++x)
    {
        sites[x] = new Site**[NY];

        for (uint y = 0; y < NY; ++y)
        {
            sites[x][y] = new Site*[NZ];

            for (uint z = 0; z < NZ; ++z)
            {
                sites[x][y][z] = new Site(x, y, z);
            }
        }
    }


    for (uint x = 0; x < NX; ++x)
    {
        for (uint y = 0; y < NY; ++y)
        {
            for (uint z = 0; z < NZ; ++z)
            {
                sites[x][y][z]->introduceNeighborhood();
            }
        }
    }

}


void KMCSolver::initializeCrystal()
{

    bool enabled = KMCDebugger_IsEnabled;
    KMCDebugger_SetEnabledTo(false);

    sites[NX/2][NY/2][NZ/2]->spawnAsFixedCrystal();
    KMCDebugger_PushTraces();

    uint crystalSizeX = round(NX*RelativeSeedSize);
    uint crystalSizeY = round(NY*RelativeSeedSize);
    uint crystalSizeZ = round(NZ*RelativeSeedSize);

    uint crystalStartX = NX/2 - crystalSizeX/2;
    uint crystalStartY = NY/2 - crystalSizeY/2;
    uint crystalStartZ = NZ/2 - crystalSizeZ/2;

    uint crystalEndX = crystalStartX + crystalSizeX;
    uint crystalEndY = crystalStartY + crystalSizeY;
    uint crystalEndZ = crystalStartZ + crystalSizeZ;

    uint solutionEndX = crystalStartX - Site::nNeighborsLimit();
    uint solutionEndY = crystalStartY - Site::nNeighborsLimit();
    uint solutionEndZ = crystalStartZ - Site::nNeighborsLimit();

    uint solutionStartX = crystalEndX + Site::nNeighborsLimit();
    uint solutionStartY = crystalEndY + Site::nNeighborsLimit();
    uint solutionStartZ = crystalEndZ + Site::nNeighborsLimit();

    for (uint i = 0; i < NX; ++i)
    {
        for (uint j = 0; j < NY; ++j)
        {
            for (uint k = 0; k < NZ; ++k)
            {

                if (i >= crystalStartX && i < crystalEndX)
                {
                    if (j >= crystalStartY && j < crystalEndY)
                    {
                        if (k >= crystalStartZ && k < crystalEndZ)
                        {
                            if (!((i == NX/2 && j == NY/2 && k == NZ/2)))
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
                    if (KMC_RNG_UNIFORM() < saturation)
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

    cout << "Initialized "
         << Site::totalActiveSites()
         << " active sites."
         << endl;

    dumpXYZ();


    KMCDebugger_SetEnabledTo(enabled);

}




void KMCSolver::getRateVariables()
{

    m_kTot = 0;
    m_accuAllRates.clear();
    m_allReactions.clear();

    Site::updateAffectedSites();

    for (uint x = 0; x < NX; ++x)
    {
        for (uint y = 0; y < NY; ++y)
        {
            for (uint z = 0; z < NZ; ++z)
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


uint KMCSolver::ptrCount = 0;
