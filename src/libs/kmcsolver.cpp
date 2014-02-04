#include "kmcsolver.h"
#include "RNG/kMCRNG.h"

#include "site.h"
#include "reactions/reaction.h"
#include "reactions/diffusion/diffusionreaction.h"

#include <sys/time.h>

#include <armadillo>
using namespace arma;

#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

KMCSolver::KMCSolver(const Setting & root)
{

    const Setting & SystemSettings = getSurfaceSetting(root, "System");
    const Setting & SolverSettings = getSurfaceSetting(root, "Solver");
    const Setting & InitializationSettings = getSurfaceSetting(root, "Initialization");

    const Setting & BoxSize = getSurfaceSetting(SystemSettings, "BoxSize");


    NX = BoxSize[0];
    NY = BoxSize[1];
    NZ = BoxSize[2];



    Site::loadNeighborLimit(SystemSettings);
    Site::setSolverPtr(this);

    Reaction::loadReactionSettings(getSurfaceSetting(root, "Reactions"));
    Reaction::setSolverPtr(this);

    DiffusionReaction::loadPotential(getSetting(root, {"Reactions", "Diffusion"}));


    nCycles = getSurfaceSetting<uint>(SolverSettings, "nCycles");

    cyclesPerOutput = getSurfaceSetting<uint>(SolverSettings, "cyclesPerOutput");


    Seed::SeedType seedType = static_cast<Seed::SeedType>(getSurfaceSetting<uint>(SolverSettings, "seedType"));

    int seed;
    switch (seedType) {
    case Seed::specific:
        seed = getSurfaceSetting<int>(SolverSettings, "specificSeed");
        break;
    case Seed::fromTime:
        seed = time(NULL);
        break;
    default:
        assert(0 == 1 && "SEED NOT SPECIFIED.");
        break;
    }

    KMC_INIT_RNG(seed);


    saturation = getSurfaceSetting<double>(InitializationSettings, "SaturationLevel");
    RelativeSeedSize = getSurfaceSetting<double>(InitializationSettings, "RelativeSeedSize");
    assert(RelativeSeedSize < 1.0 && "The seed size cannot exceed the box size.");


    sites = new Site***[NX];

    for (uint x = 0; x < NX; ++x) {

        sites[x] = new Site**[NY];
        for (uint y = 0; y < NY; ++y) {

            sites[x][y] = new Site*[NZ];

            for (uint z = 0; z < NZ; ++z) {

                sites[x][y][z] = new Site(x, y, z);
            }
        }
    }


    for (uint x = 0; x < NX; ++x) {
        for (uint y = 0; y < NY; ++y) {
            for (uint z = 0; z < NZ; ++z) {
                sites[x][y][z]->introduceNeighborhood();
            }
        }
    }



}

KMCSolver::~KMCSolver()
{
    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {

                delete sites[i][j][k];
            }

            delete [] sites[i][j];
        }

        delete [] sites[i];
    }

    delete [] sites;

    Site::resetAll();

}



void KMCSolver::run(){

    uint choice;
    double R;

    initialize();

    while(cycle < nCycles)
    {

        getRateVariables();

        R = kTot*KMC_RNG_UNIFORM();

        choice = getReactionChoice(R);

        allReactions[choice]->execute();

        if (cycle%cyclesPerOutput == 0)
        {
            dumpOutput();
            dumpXYZ();
        }


        totalTime += Reaction::getScale()/kTot;
        cycle++;

    }

}



void KMCSolver::dumpXYZ()
{
    stringstream s;
    s << "kMC" << outputCounter++ << ".xyz";

    ofstream o;
    o.open("outfiles/" + s.str());
    o << Site::totalActiveSites() << "\n - ";

    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {
                if (sites[i][j][k]->active()) {
                    o << "\nC " << i << " " << j << " " << k << " " << sites[i][j][k]->nNeighbors();
                }
            }
        }
    }

    o.close();

}

void KMCSolver::dumpOutput()
{
    cout << setw(5) << right << setprecision(1) << fixed
         << (double)cycle/nCycles*100 << "%   "
         << outputCounter
         << endl;
}

void KMCSolver::getAllNeighbors()
{
    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {
                sites[i][j][k]->countNeighbors();
            }
        }
    }
}


void KMCSolver::setDiffusionReactions()
{

    Site* currentSite;
    Site* destination;

    //Loop over all sites
    for (uint x = 0; x < NX; ++x) {
        for (uint y = 0; y < NY; ++y) {
            for (uint z = 0; z < NZ; ++z) {

                currentSite = sites[x][y][z];

                //For each site, loop over all neightbours
                for (uint i = 0; i < 3; ++i) {
                    for (uint j = 0; j < 3; ++j) {
                        for (uint k = 0; k < 3; ++k) {

                            destination = currentSite->getNeighborhood()[Site::nNeighborsLimit() + i - 1][Site::nNeighborsLimit() + j - 1][Site::nNeighborsLimit() + k - 1];

                            //This menas we are at the current site.
                            if(destination == currentSite) {
                                assert((i == 1) && (j == 1) && (k == 1));
                                continue;
                            }

                            //And add diffusion reactions
                            currentSite->addReaction(new DiffusionReaction(destination));

                        }
                    }
                }

                //Then we update the site reactions based on the current setup
                currentSite->updateReactions();

                //And calculate the process rates (only dependent on other sites, not reactions)
                currentSite->calculateRates();

            }
        }
    }
}


void KMCSolver::initialize()
{

    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {
                if (KMC_RNG_UNIFORM() > 1 - saturation) {

                    sites[i][j][k]->activate();

                }
            }
        }
    }


    int D = NX*RelativeSeedSize/2;

    for (uint i = NX/2 - D; i < NX/2 + D; ++i) {
        for (uint j = NY/2 - D; j < NY/2 + D; ++j) {
            for (uint k = NZ/2 - D; k < NZ/2 + D; ++k) {
                if (!sites[i][j][k]->active()){
                    sites[i][j][k]->activate();
                }
            }
        }
    }

//    double E0 = Site::totalEnergy();
//    ucube alln(NX, NY, NZ);
//    ucube alln2(NX, NY, NZ);

//    for (uint i = 0; i < NX; ++i) {
//        for (uint j = 0; j < NY; ++j) {
//            for (uint k = 0; k < NZ; ++k) {
//                alln(i, j, k) = sites[i][j][k]->nNeighbors();
//                alln2(i, j, k) = sites[i][j][k]->nNeighbors(1);

//                sites[i][j][k]->reset();


//            }
//        }
//    }

//    double isZero = Site::totalEnergy();

//    getAllNeighbors();

//    for (uint i = 0; i < NX; ++i) {
//        for (uint j = 0; j < NY; ++j) {
//            for (uint k = 0; k < NZ; ++k) {
//                assert(alln(i, j, k) == sites[i][j][k]->nNeighbors());
//                assert(alln2(i, j, k) == sites[i][j][k]->nNeighbors(1));


//            }
//        }
//    }

//    double E1 = Site::totalEnergy();

//    cout << E0 << "  " << E1 << "  " << isZero << endl;
//    assert(fabs(E0 - E1) < 1E-3);

    cout << "Initialized "
         << Site::totalActiveSites()
         << " active sites."
         << endl;

    dumpXYZ();

    setDiffusionReactions();

}



void KMCSolver::getRateVariables()
{

    kTot = 0;
    accuAllRates.clear();
    allReactions.clear();

    for (uint x = 0; x < NX; ++x) {
        for (uint y = 0; y < NY; ++y) {
            for (uint z = 0; z < NZ; ++z) {
                for (Reaction* reaction : sites[x][y][z]->activeReactions()) {
                    kTot += reaction->rate();
                    accuAllRates.push_back(kTot);
                    allReactions.push_back(reaction);
                }
            }
        }
    }

}


uint KMCSolver::getReactionChoice(double R)
{

    uint imax = accuAllRates.size() - 1;
    uint MAX = imax;
    uint imin = 0;
    uint imid = 1;

    // continue searching while imax != imin + 1
    while (imid != imin)
    {
        if (imin > imax) {
            cout << "caught: " << imin << " " << imid << " " << imax << endl;
            return imin;
        }

        // calculate the midpoint for roughly equal partition
        imid = imin + (imax - imin)/2;

        //Is the upper limit above mid?
        if (R > accuAllRates.at(imid))
        {

            if (imid == MAX) {
                return MAX;
            }

            //Are we just infront of the limit?
            else if (R < accuAllRates.at(imid + 1))
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

    cout << "binary: " << imid << "  " << imin << "  " << imax << endl;
//    assert(imax == imid + 1 && "LOCATING RATE FAILED");
    return imid + 1;

}


