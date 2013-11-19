#include "kmcsolver.h"
#include "RNG/kMCRNG.h"

#include "site.h"
#include "reactions/reaction.h"

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

                sites[x][y][z] = new Site(x, y, z);
            }
        }
    }

}

void KMCSolver::addReaction(Reaction *reaction, uint &x, uint &y, uint &z)
{
    reaction->setKMCSolverObject(this);
    sites[x][y][z]->addReaction(reaction);
}

void KMCSolver::run(){

    nTot = 0;
    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {
                if ((i < NX/10) || (i > 9*NX/10)) {
                    sites[i][j][k]->activate();
                    nTot++;
                } else {
                    if (KMC_RNG_UNIFORM() > 0.85) {

                        sites[i][j][k]->activate();
                        nTot++;
                    }
                }
            }
        }
    }

    dumpXYZ();

    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {
                getNeighbours(i, j, k);
            }
        }
    }

    double T = 1E7;

    double dt;

    std::vector<double> allRates;
    std::vector<std::vector<uint>> transitions;

    double kTot;
    double Enn = 0.05;
    double Ennn = 0.01;
    double temperature = 300;
    double mu = 1;
    while(counter < 100000) {

        allRates.clear();
        kTot = 0;
        transitions.clear();
        for (uint i = 0; i < NX; ++i) {
            for (uint j = 0; j < NY; ++j) {
                for (uint k = 0; k < NZ; ++k) {

                    if(!sites[i][j][k]->active()) {
                        continue;
                    }

                    uint nn = neighbours(i, j)(k).n_rows;
                    uint nnn = nextNeighbours(i, j)(k).n_rows;

                    double Eijk = nn*Enn + nnn*Ennn;

                    for (uint ci = 0; ci < 3; ++ci) {
                        uint inp = (i + delta(ci) + NX)%NX;

                        for (uint cj = 0; cj < 3; ++cj) {
                            uint jnp = (j + delta(cj) + NY)%NY;

                            for (uint ck = 0; ck < 3; ++ck) {

                                if((ci == 1) && (cj == 1) && (ck == 1)) {
                                    continue;
                                }
                                uint knp = (k + delta(ck)+ NZ)%NZ;


                                if (!sites[inp][jnp][knp]->active()) {
                                    uint nnsp = neighbours(inp, jnp)(knp).n_rows;
                                    uint nnnsp = nextNeighbours(inp, jnp)(knp).n_rows;
                                    double Esp = 1.5*((nn + nnsp)*Enn + (nnn + nnnsp)*Ennn);

                                    double rate = mu*exp(-(Eijk-Esp)/temperature);
                                    kTot += rate;
                                    allRates.push_back(kTot);
                                    std::vector<uint> trans = {i, j, k, inp, jnp, knp};
                                    transitions.push_back(trans);
                                }
                            }
                        }
                    }

                }
            }
        }

        double R = kTot*KMC_RNG_UNIFORM();

        int choice = 0;
        double r = 0;
        while(allRates.at(choice) < R) {
            choice++;
        }

        std::vector<uint> chosen = transitions.at(choice);
        uint x0 = chosen.at(0);
        uint y0 = chosen.at(1);
        uint z0 = chosen.at(2);
        uint x1 = chosen.at(3);
        uint y1 = chosen.at(4);
        uint z1 = chosen.at(5);

        deactivateSite(x0,y0, z0);
        activateSite(x1, y1, z1);

        dt = 1.0/kTot;
        t += dt;

        counter2++;

        if (counter2%200 == 0){
            dumpXYZ();
            counter++;
            cout << t/T << endl;
        }

    }

}


void KMCSolver::test()
{


    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {
                if (KMC_RNG_UNIFORM() > 0.75) {

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


    cout << nTot/(double) (NX*NY*NZ) << endl;

    double T = 1E7;

    double dt;

    uint i, j, k;

    dt = 1;

    uint nV, n;
    double E0, E0i;
    while(t < T) {

        i = (uint) (NX*KMC_RNG_UNIFORM());
        j = (uint) (NY*KMC_RNG_UNIFORM());
        k = (uint) (NZ*KMC_RNG_UNIFORM());

        n = neighbours(i, j)(k).n_rows;
        nV = vacantNeighbours(i, j)(k).n_rows;

        E0 = pow(n/26.0, 1.015);

        if (sites[i][j][k]->active()) {
            double rateDiff = 1-E0;
            double rateAbsorbed = E0;
            double totalRate = rateDiff + rateAbsorbed;

            double R = KMC_RNG_UNIFORM()*totalRate;

            if (R < E0) {
                reactionDiffusion(i, j, k);
            } else {
                deactivateSite(i, j, k);
            }

        } else {
            double rateReleased = 1-E0;
            double ratePass = E0;
            double totalRate = rateReleased + ratePass;

            double R = KMC_RNG_UNIFORM()*totalRate;

            if (R < E0) {
                activateSite(i, j, k);
            }
        }

        //        for (uint x = 0; x < NX; ++x) {
        //            for (uint y = 0; y < NY; ++y) {
        //                for (uint z = 0; z < NZ; ++z) {
        ////                    getNeighbours(x, y, z);
        //                    umat & A = vacantNeighbours(x, y)(z);
        //                    umat & B = neighbours(x, y)(z);

        //                    int X = vacantNeighbours(x, y)(z).n_elem/3 + neighbours(x, y)(z).n_elem/3;
        //                    if (X != 26) {
        //                        cout << "FAIL AT UPDATING NEIGBOURS " << X << endl;
        //                        cout << A << endl << B << endl;
        //                        exit(1);
        //                    }
        //                    for (int s = 0; s < vacantNeighbours(x, y)(z).n_rows; ++s) {
        //                        if ((x == A(s, 0)) && (y == A(s, 1)) && (z == A(s, 2))) {
        //                            cout << "FAIL. CONTAINED SELF" << endl;
        //                            exit(1);
        //                        }
        //                    }
        //                }
        //            }
        //        }

        t += dt;

        counter2++;

        if (counter2%1000 == 0){
            dumpXYZ();
            counter++;
            cout << t/T << endl;
        }

    }
}



void KMCSolver::dumpXYZ()
{
    ofstream o;
    stringstream s;
    s << "kMC" << counter << ".xyz";
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

void KMCSolver::updateNextNeighbour(uint & x, uint& y, uint &z, const urowvec & newRow, bool activate)
{

    //activate=true means the particle at (i j k) needs to be added to all surrounding neighbour lists
    if (activate) {
        nextNeighbours(x, y)(z).insert_rows(0, newRow);
    } else {
        for (uint l = 0; l < nextNeighbours(x, y)(z).n_rows; ++l) {

            uint & i = nextNeighbours(x, y)(z)(l, 0);
            uint & j = nextNeighbours(x, y)(z)(l, 1);
            uint & k = nextNeighbours(x, y)(z)(l, 2);

            if ((x == i) && (y==j) && (z == k)) {
                nextNeighbours(x, y)(z).shed_row(l);
                return;
            }
        }
    }
}

void KMCSolver::reactionDiffusion(uint i, uint j, uint k)
{

    uint nVacant = vacantNeighbours(i, j)(k).n_rows;
    uint loc     = (uint)(nVacant*KMC_RNG_UNIFORM());

    uint ip = vacantNeighbours(i, j)(k)(loc, 0);
    uint jp = vacantNeighbours(i, j)(k)(loc, 1);
    uint kp = vacantNeighbours(i, j)(k)(loc, 2);

    activateSite(ip, jp, kp);
    deactivateSite(i, j, k);

}

void KMCSolver::activateSite(uint i, uint j, uint k)
{

    //create the new site
    sites[i][j][k]->activate();
    nTot++;

    updateNeighbourLists(neighbours, vacantNeighbours, i, j, k);

}

void KMCSolver::deactivateSite(uint i, uint j, uint k)
{

    sites[i][j][k]->deactivate();
    nTot--;

    updateNeighbourLists(vacantNeighbours, neighbours, i, j, k);

}

void KMCSolver::updateNeighbourLists(field<field<umat>> & A, field<field<umat>> & B,
                                     uint i, uint j, uint k, bool activate)
{
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

                //                if ((abs((int)(k - knp))%NZ != 1) && (abs((int)(k - knp))%NZ != 29) &&(abs((int)(j - jnp))%NY != 1)&& (abs((int)(j - jnp))%NY != 29)  &&(abs((int)(i - inp))%NX != 1) && (abs((int)(i - inp))%NX != 29) ){
                //                    cout << "MOVED MORE THAN ONE ";
                //                    cout << i << "  " << inp << endl;
                //                    cout << j << "  " << jnp << endl;
                //                    cout << k << "  " << knp << endl;
                //                    exit(1);
                //                }

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

    int C = 0;

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
            C+=2;

        }
    }



    //LEFT AND RIGHT LAYER
    for (uint ci = 0; ci < 5; ++ci){

        X = (X0 + ci + NX)%NX;
        for (uint ck = 1; ck < 4; ++ck) {

            Z = (Z0 + ck + NZ)%NZ;

            updateNextNeighbour(X, Y0, Z, newRow, activate);
            updateNextNeighbour(X, Y1, Z, newRow, activate);

            C+=2;
        }
    }


    //BACK AND FRONT LAYER
    for (uint cj = 1; cj < 4; ++cj){

        Y = (Y0 + cj + NY)%NY;
        for (uint ck = 1; ck < 4; ++ck) {

            Z = (Z0 + ck + NZ)%NZ;

            updateNextNeighbour(X0, Y, Z, newRow, activate);
            updateNextNeighbour(X1, Y, Z, newRow, activate);

            C+=2;
        }
    }


//    cout << C << endl;

}
