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

    neighbours.set_size(NX, NY);
    nextNeighbours.set_size(NX, NY);
    vacantNeighbours.set_size(NX, NY);

    sites = new Site***[NX];

    for (uint x = 0; x < NX; ++x) {

        sites[x] = new Site**[NY];
        for (uint y = 0; y < NY; ++y) {

            neighbours(x, y).set_size(NZ);
            nextNeighbours(x, y).set_size(NZ);
            vacantNeighbours(x, y).set_size(NZ);
            sites[x][y] = new Site*[NZ];

            for (uint z = 0; z < NZ; ++z) {

                sites[x][y][z] = new Site(x, y, z, this);
            }
        }
    }

}

void KMCSolver::run(){

    uint choice;

    nTot = 0;

    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {
                if (KMC_RNG_UNIFORM() > 0.9) {

                    sites[i][j][k]->activate();
                    nTot++;
                }
            }
        }
    }


    int D = 10;
    for (uint i = NX/2 - D; i < NX/2 + D; ++i) {
        for (uint j = NY/2 - D; j < NY/2 + D; ++j) {
            for (uint k = NZ/2 - D; k < NZ/2 + D; ++k) {
                if (!sites[i][j][k]->active()){
                    sites[i][j][k]->activate();
                    nTot++;
                }
            }
        }
    }

    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {
                getNeighbours(i, j, k);
            }
        }
    }

    cout << nTot << endl;

    dumpXYZ();

    setDiffusionReactions();

//    wall_clock wc;
//    double t1 = 0;
    while(counter < 100000) {
//        wc.tic();
        getRateVariables();

        double R = kTot*KMC_RNG_UNIFORM();

        choice = getReactionChoice(R);

        Reaction* chosenReaction = allReactions[choice];

        reactionAffectedSites.clear();
        chosenReaction->execute();

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

    o << nTot << "\n - ";
    uint COUNT = 0;
    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {
                if (sites[i][j][k]->active()) {
                    o << "\nC " << i << " " << j << " " << k << " " << neighbours(i, j)(k).n_rows;
                    COUNT++;
                }
            }
        }
    }
    if (COUNT != nTot) {
        cout << "FAIL FAIL FAIL "<< COUNT << "  " << nTot << endl;
        exit(1);
    }
    o.close();

}

void KMCSolver::getNeighbours(uint i, uint j, uint k)
{
    //Count neighbours
    uint n = 0;

    //Count vacant neighbours
    uint nf = 0;

    //Count next neighbours
    uint nn = 0;

    //Create maximal size matrices. Reshape them based on above integers post.
    umat localNeighbours(27, 3);
    umat localVacantNeighbours(27, 3);
    umat localNextNeighbours(98, 3);



    //Find closest neighbours and vacancies
    for (uint ci = 0; ci < 3; ++ci) {
        uint inp = (i + delta(ci) + NX)%NX;

        for (uint cj = 0; cj < 3; ++cj) {
            uint jnp = (j + delta(cj) + NY)%NY;

            for (uint ck = 0; ck < 3; ++ck) {

                if((ci == 1) && (cj == 1) && (ck == 1)) {
                    continue;
                }

                uint knp = (k + delta(ck) + NZ)%NZ;


                if (sites[inp][jnp][knp]->active()) {
                    localNeighbours(n, 0) = inp;
                    localNeighbours(n, 1) = jnp;
                    localNeighbours(n, 2) = knp;
                    n++;
                } else {
                    localVacantNeighbours(nf, 0) = inp;
                    localVacantNeighbours(nf, 1) = jnp;
                    localVacantNeighbours(nf, 2) = knp;
                    nf++;
                }

            }
        }
    }




    //Find next neighbours
    uint X, Y, Z;
    uint X0 = (i - 2 + NX) % NX;
    uint Y0 = (j - 2 + NY) % NY;
    uint Z0 = (k - 2 + NZ) % NZ;

    uint X1 = (i + 2 + NX) % NX;
    uint Y1 = (j + 2 + NY) % NY;
    uint Z1 = (k + 2 + NZ) % NZ;

    //TOP AND BOTTOM LAYER
    for (uint ci = 0; ci < 5; ++ci) {

        X = (X0 + ci + NX)%NX;
        for (uint cj = 0; cj < 5; ++cj) {

            Y = (Y0 + cj + NY)%NY;

            if (sites[X][Y][Z0]->active()){
                localNextNeighbours(nn, 0) = X;
                localNextNeighbours(nn, 1) = Y;
                localNextNeighbours(nn, 2) = Z0;

                nn++;
            }

            if (sites[X][Y][Z1]->active()){
                localNextNeighbours(nn, 0) = X;
                localNextNeighbours(nn, 1) = Y;
                localNextNeighbours(nn, 2) = Z1;

                nn++;
            }

        }
    }



    //LEFT AND RIGHT LAYER
    for (uint ci = 0; ci < 5; ++ci){

        X = (X0 + ci + NX)%NX;
        for (uint ck = 1; ck < 4; ++ck) {

            Z = (Z0 + ck + NZ)%NZ;

            if (sites[X][Y0][Z]->active()){
                localNextNeighbours(nn, 0) = X;
                localNextNeighbours(nn, 1) = Y0;
                localNextNeighbours(nn, 2) = Z;

                nn++;
            }

            if (sites[X][Y1][Z]->active()){
                localNextNeighbours(nn, 0) = X;
                localNextNeighbours(nn, 1) = Y1;
                localNextNeighbours(nn, 2) = Z;

                nn++;
            }
        }
    }


    //BACK AND FRONT LAYER
    for (uint cj = 1; cj < 4; ++cj){

        Y = (Y0 + cj + NY)%NY;
        for (uint ck = 1; ck < 4; ++ck) {

            Z = (Z0 + ck + NZ)%NZ;

            if (sites[X0][Y][Z]->active()){
                localNextNeighbours(nn, 0) = X0;
                localNextNeighbours(nn, 1) = Y;
                localNextNeighbours(nn, 2) = Z;

                nn++;
            }

            if (sites[X1][Y][Z]->active()){
                localNextNeighbours(nn, 0) = X1;
                localNextNeighbours(nn, 1) = Y;
                localNextNeighbours(nn, 2) = Z;

                nn++;
            }

        }
    }



    //Resize matrices mased on the integer count value and fill the neighbour list.
    if (n != 0) {
        neighbours(i, j)(k) = localNeighbours(span(0, n-1), span::all);
    } else {
        neighbours(i, j)(k).reset();
    }

    if (nf != 0) {
        vacantNeighbours(i, j)(k) = localVacantNeighbours(span(0, nf-1), span::all);
    } else {
        vacantNeighbours(i, j)(k).reset();
    }

    if (nn != 0) {
        nextNeighbours(i, j)(k) = localNextNeighbours(span(0, nn-1), span::all);
    } else {
        nextNeighbours(i, j)(k).reset();
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

void KMCSolver::updateNextNeighbour(uint & x, uint& y, uint &z, const urowvec & newRow, bool activate)
{


    bool notUpdated = pushToRateQueue(sites[x][y][z]);

    if (!notUpdated) {
        return;
    }

    uint xnew = newRow(0);
    uint ynew = newRow(1);
    uint znew = newRow(2);


    //activate=true means the particle at (i j k) needs to be added to all surrounding neighbour lists
    if (activate) {
        nextNeighbours(x, y)(z).insert_rows(0, newRow);
    } else {
        for (uint l = 0; l < nextNeighbours(x, y)(z).n_rows; ++l) {

            uint & i = nextNeighbours(x, y)(z)(l, 0);
            uint & j = nextNeighbours(x, y)(z)(l, 1);
            uint & k = nextNeighbours(x, y)(z)(l, 2);

            if ((xnew == i) && (ynew==j) && (znew == k)) {
                nextNeighbours(x, y)(z).shed_row(l);
                return;
            }
        }
    }
}

void KMCSolver::activateSite(Site* site)
{

    if (site->active()) {
        cout << "fail. activating active site... rate should be zero." << endl;
        exit(1);
    }


    //create the new site
    site->activate();
    nTot++;

    updateNeighbourLists(neighbours, vacantNeighbours, site, true);

}

void KMCSolver::deactivateSite(Site* site)
{

    if (!site->active()) {
        cout << "fail2. Deactivating deactivated site.. rate should be zero." << endl;
        exit(1);
    }

    site->deactivate();
    nTot--;

    updateNeighbourLists(vacantNeighbours, neighbours, site);

}

void KMCSolver::updateNeighbourLists(field<field<umat>> & A, field<field<umat>> & B,
                                     Site* site, bool activate)
{

    uint i = site->x();
    uint j = site->y();
    uint k = site->z();

    urowvec newRow = {i, j, k};

    //All comments works mirrored if we remove a particle

    //loop over the box surrounding the new site and
    //add the site as neighbours to surrounding sites
    //and then also remove it from the vacancy lists
    for (uint ci = 0; ci < 3; ++ci) {
        uint inp = (i + delta(ci) + NX)%NX;

        for (uint cj = 0; cj < 3; ++cj) {
            uint jnp = (j + delta(cj) + NY)%NY;

            for (uint ck = 0; ck < 3; ++ck) {

                if((ci == 1) && (cj == 1) && (ck == 1)) {
                    continue;
                }

                uint knp = (k + delta(ck) + NZ)%NZ;


                bool notUpdated = pushToRateQueue(sites[inp][jnp][knp]);

                if (!notUpdated) {
                    continue;
                }

                //add the new site as neighbour to surrounding site
                A(inp, jnp)(knp).insert_rows(0, newRow);

                //remove the site from the vacancy list of the surrounding site
                for(uint l = 0; l < B(inp, jnp)(knp).n_rows; ++l){

                    uint & x = B(inp, jnp)(knp)(l, 0);
                    uint & y = B(inp, jnp)(knp)(l, 1);
                    uint & z = B(inp, jnp)(knp)(l, 2);

                    if ((x == i) && (y==j) && (z == k)) {
                        B(inp, jnp)(knp).shed_row(l);
                        break;
                    }
                    //                    if (l == B(inp, jnp)(knp).n_rows - 1) {
                    //                        cout << "SHOULD NEVER GET HERE" << endl;
                    //                        cout << i << " " << j << " " << k << endl;
                    //                        cout <<     B(inp, jnp)(knp) << endl;
                    //                        exit(1);
                    //                    }
                }


            }
        }
    }


    uint X, Y, Z;
    uint X0 = (i - 2 + NX) % NX;
    uint Y0 = (j - 2 + NY) % NY;
    uint Z0 = (k - 2 + NZ) % NZ;

    uint X1 = (i + 2 + NX) % NX;
    uint Y1 = (j + 2 + NY) % NY;
    uint Z1 = (k + 2 + NZ) % NZ;

    //TOP AND BOTTOM LAYER
    for (uint ci = 0; ci < 5; ++ci) {

        X = (X0 + ci + NX)%NX;
        for (uint cj = 0; cj < 5; ++cj) {

            Y = (Y0 + cj + NY)%NY;

            updateNextNeighbour(X, Y, Z0, newRow, activate);
            updateNextNeighbour(X, Y, Z1, newRow, activate);

        }
    }



    //LEFT AND RIGHT LAYER
    for (uint ci = 0; ci < 5; ++ci){

        X = (X0 + ci + NX)%NX;
        for (uint ck = 1; ck < 4; ++ck) {

            Z = (Z0 + ck + NZ)%NZ;

            updateNextNeighbour(X, Y0, Z, newRow, activate);
            updateNextNeighbour(X, Y1, Z, newRow, activate);

        }
    }


    //BACK AND FRONT LAYER
    for (uint cj = 1; cj < 4; ++cj){

        Y = (Y0 + cj + NY)%NY;
        for (uint ck = 1; ck < 4; ++ck) {

            Z = (Z0 + ck + NZ)%NZ;

            updateNextNeighbour(X0, Y, Z, newRow, activate);
            updateNextNeighbour(X1, Y, Z, newRow, activate);

        }
    }


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

