#include "kmcsolver.h"
#include "RNG/kMCRNG.h"

#include "site.h"
#include "reactions/reaction.h"
#include "reactions/diffusion/diffusionreaction.h"

#include <sys/time.h>

#include <armadillo>
using namespace arma;

#include <iostream>
#include <fstream>

using namespace std;

KMCSolver::KMCSolver(uint NX, uint NY, uint NZ) :
    NX(NX),
    NY(NY),
    NZ(NZ)
{

    KMC_INIT_RNG(time(NULL));

    for (uint i = 0; i < Site::neighborhoodLength; ++i) {
        for (uint j = 0; j < Site::neighborhoodLength; ++j) {
            for (uint k = 0; k < Site::neighborhoodLength; ++k) {

                Site::levelMatrix(i, j, k) = Site::getLevel(std::abs((int)i - (int)Site::nNeighborsLimit),
                                                            std::abs((int)j - (int)Site::nNeighborsLimit),
                                                            std::abs((int)k - (int)Site::nNeighborsLimit));
            }
        }
    }

    sites = new Site***[NX];

    for (uint x = 0; x < NX; ++x) {

        sites[x] = new Site**[NY];
        for (uint y = 0; y < NY; ++y) {

            sites[x][y] = new Site*[NZ];

            for (uint z = 0; z < NZ; ++z) {

                sites[x][y][z] = new Site(x, y, z, this);
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



void KMCSolver::run(){

    uint choice;

    initialize();

//    wall_clock wc;
//    double t1 = 0;
    while(counter < 100000) {
//        wc.tic();

        getAllNeighbors();
        getRateVariables();

        double R = kTot*KMC_RNG_UNIFORM();

        choice = getReactionChoice(R);

        Reaction* chosenReaction = allReactions[choice];

        reactionAffectedSites.clear();
        chosenReaction->execute();



        getAllNeighbors();

        updateRates();



        counter2++;

        if (counter2%250 == 0){
            cout << counter << endl;
            dumpXYZ();
        }


        t += 1.0/kTot;
//        t1 += wc.toc();

    }

}



void KMCSolver::dumpXYZ()
{
    ofstream o;
    stringstream s;
    s << "kMC" << counter++ << ".xyz";
    o.open("outfiles/" + s.str());

    o << Site::totalActiveSites << "\n - ";
    uint COUNT = 0;
    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {
                if (sites[i][j][k]->active()) {
                    o << "\nC " << i << " " << j << " " << k << " " << sites[i][j][k]->nNeighbors();
                    COUNT++;
                }
            }
        }
    }
    if (COUNT != Site::totalActiveSites) {
        cout << "FAIL FAIL FAIL "<< COUNT << "  " << Site::totalActiveSites << endl;
        exit(1);
    }
    o.close();

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
    //Loop over all sites
    for (uint x = 0; x < NX; ++x) {
        for (uint y = 0; y < NY; ++y) {
            for (uint z = 0; z < NZ; ++z) {

                //For each site, loop over all neightbours
                for (uint dx_i = 0; dx_i < 3; ++dx_i) {
                    uint x1 = (x + delta(dx_i) + NX)%NX;

                    for (uint dy_i = 0; dy_i < 3; ++dy_i) {
                        uint y1 = (y + delta(dy_i) + NY)%NY;

                        for (uint dz_i = 0; dz_i < 3; ++dz_i) {

                            //This menas we are at the current site.
                            if((dx_i == 1) && (dy_i == 1) && (dz_i == 1)) {
                                continue;
                            }
                            uint z1 = (z + delta(dz_i)+ NZ)%NZ;

                            //And add diffusion reactions
                            sites[x][y][z]->addReaction(new DiffusionReaction(sites[x1][y1][z1]));

                        }
                    }
                }

                //Then we update the site reactions based on the current setup
                sites[x][y][z]->updateReactions();

                //And calculate the process rates (only dependent on other sites, not reactions)
                sites[x][y][z]->calculateRates();

            }
        }
    }
}


void KMCSolver::initialize()
{

    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {
                if (KMC_RNG_UNIFORM() > 0.98) {

                    sites[i][j][k]->activate();

                }
            }
        }
    }


    int D = NX/5;
    for (uint i = NX/2 - D; i < NX/2 + D; ++i) {
        for (uint j = NY/2 - D; j < NY/2 + D; ++j) {
            for (uint k = NZ/2 - D; k < NZ/2 + D; ++k) {
                if (!sites[i][j][k]->active()){
                    sites[i][j][k]->activate();
                }
            }
        }
    }


    getAllNeighbors();

    cout << Site::totalActiveSites << endl;

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


void KMCSolver::updateRates()
{

    for (Site* affectedSite : reactionAffectedSites) {
        affectedSite->updateReactions();
        affectedSite->calculateRates();
    }

}

Reaction* KMCSolver::getChosenReaction(uint choice)
{
    uint K = 0;

    for (uint x = 0; x < NX; ++x) {
        for (uint y = 0; y < NY; ++y) {
            for (uint z = 0; z < NZ; ++z) {

                for (Reaction* reaction : sites[x][y][z]->activeReactions()) {

                    if (K == choice) {
                        return reaction;
                    }

                    K++;
                }

            }
        }
    }

    cout << "FAIL AT CHOOSING REACTION" << endl;
    exit(1);
}

uint KMCSolver::getReactionChoice(double R)
{

    uint imax = accuAllRates.size();
    uint imin = 0;
    uint imid = 1;

    // continue searching while imax != imin + 1
    while (imid != imin)
    {
        // calculate the midpoint for roughly equal partition
        imid = imin + (imax - imin)/2;

        if (imid == imin) {
            return imid;
        }

        //perform two tests due to integer division
        else if (accuAllRates.at(imid - 1) < R && accuAllRates.at(imid) > R) {
            return imid-1;
        }
        else if (accuAllRates.at(imid) < R && accuAllRates.at(imid + 1) > R) {
            return imid;
        }

        else if (accuAllRates.at(imid) < R) {
            imin = imid + 1;
        } else {
            imax = imid - 1;
        }


    }

    cout << "FAILED" << endl;
    return 0;

}

bool KMCSolver::pushToRateQueue(Site *affectedSite)
{
    for (const Site * prevAffectedSite : reactionAffectedSites) {

        //Skip out of the function in the site is already queued
        if (prevAffectedSite == affectedSite) {
            return false;
        }
    }

    reactionAffectedSites.push_back(affectedSite);

    return true;
}

