#include "solver.h"
#include "RNG/kMCRNG.h"

#include <sys/time.h>

#include <armadillo>
using namespace arma;

#include <iostream>
#include <fstream>

using namespace std;

BoolCube makeBoolCube(uint N, uint M, uint L, bool fill=false) {
    BoolCube boolCube = new bool**[N];

    for (uint i = 0; i < N; ++i) {
        boolCube[i] = new bool*[M];

        for (uint j = 0; j < M; ++j) {
            boolCube[i][j] = new bool[L];

            for (uint k = 0; k < L; ++k) {
                boolCube[i][j][k] = fill;
            }

        }
    }

    return boolCube;
}

Solver::Solver()
{

    KMC_INIT_RNG(time(NULL));

    siteCube = makeBoolCube(N, M, L);
    neighbours.set_size(N, M);
    vacantNeighbours.set_size(N, M);

    for (uint i = 0; i < N; ++i) {
        for (uint j = 0; j < M; ++j) {
            neighbours(i, j).set_size(L);
            vacantNeighbours(i, j).set_size(L);

            for (uint k = 0; k < L; ++k) {
                if (KMC_RNG_UNIFORM() > 0.75) {

                    siteCube[i][j][k] = true;
                    nTot++;
                }
            }
        }
    }

    for (uint i = 0; i < N; ++i) {
        for (uint j = 0; j < M; ++j) {
            for (uint k = 0; k < L; ++k) {
                getNeighbours(i, j, k);
            }
        }
    }


    cout << nTot/(double) (N*M*L) << endl;

    double T = 1E7;

    double dt;

    uint i, j, k;

    dt = 1;

    uint nV, n;
    double E0, E0i;
    while(t < T) {

        i = (uint) (N*KMC_RNG_UNIFORM());
        j = (uint) (M*KMC_RNG_UNIFORM());
        k = (uint) (L*KMC_RNG_UNIFORM());

        n = neighbours(i, j)(k).n_rows;
        nV = vacantNeighbours(i, j)(k).n_rows;

        E0 = getRate(pow(n/26.0, 1 + beta));

        if (siteCube[i][j][k]) {
            double rateDiff = 1-E0;
            double rateAbsorbed = E0;
            double totalRate = rateDiff + rateAbsorbed;

            double R = KMC_RNG_UNIFORM()*totalRate;

            if (R < E0) {
                reactionDiffusion(i, j, k);
            } else {
                reactionDeletion(i, j, k);
            }

        } else {
            double rateReleased = 1-E0;
            double ratePass = E0;
            double totalRate = rateReleased + ratePass;

            double R = KMC_RNG_UNIFORM()*totalRate;

            if (R < E0) {
                reactionCreation(i, j, k);
            }
        }

        //        for (uint x = 0; x < N; ++x) {
        //            for (uint y = 0; y < M; ++y) {
        //                for (uint z = 0; z < L; ++z) {
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
            dump();
            counter++;
            cout << t/T << endl;
        }

    }
}

void Solver::test()
{
    cout << "test" << endl;
}

void Solver::dump()
{
    ofstream o;
    stringstream s;
    s << "kMC" << counter << ".xyz";
    o.open("outfiles/" + s.str());

    o << nTot << "\n - ";
    for (uint i = 0; i < N; ++i) {
        for (uint j = 0; j < M; ++j) {
            for (uint k = 0; k < L; ++k) {
                if (siteCube[i][j][k]) {
                    o << "\nC " << i << " " << j << " " << k << " " << neighbours(i, j)(k).n_rows;
                }
            }
        }
    }
    o.close();

}

double Solver::getRate(double E)
{
    return E;
}

void Solver::getNeighbours(uint i, uint j, uint k)
{
    //Count neighbours
    uint n = 0;

    //Count vacant neighbours
    uint nf = 0;

    umat localNeighbours(27, 3);
    umat localVacantNeighbours(27, 3);
    for (uint ci = 0; ci < 3; ++ci) {
        uint inp = (i + delta(ci) + N)%N;

        for (uint cj = 0; cj < 3; ++cj) {
            uint jnp = (j + delta(cj) + M)%M;

            for (uint ck = 0; ck < 3; ++ck) {

                if((ci == 1) && (cj == 1) && (ck == 1)) {
                    continue;
                }

                uint knp = (k + delta(ck) + L)%L;


                if (siteCube[inp][jnp][knp]) {
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

}

void Solver::reactionDiffusion(uint i, uint j, uint k)
{

    uint nVacant = vacantNeighbours(i, j)(k).n_rows;
    uint loc     = (uint)(nVacant*KMC_RNG_UNIFORM());

    uint ip = vacantNeighbours(i, j)(k)(loc, 0);
    uint jp = vacantNeighbours(i, j)(k)(loc, 1);
    uint kp = vacantNeighbours(i, j)(k)(loc, 2);

    reactionCreation(ip, jp, kp);
    reactionDeletion(i, j, k);

}

void Solver::reactionCreation(uint i, uint j, uint k)
{

    //create the new site
    siteCube[i][j][k] = true;
    nTot++;

    updateNeighbourLists(neighbours, vacantNeighbours, i, j, k);

}

void Solver::reactionDeletion(uint i, uint j, uint k)
{
    //delete the site
    siteCube[i][j][k] = false;
    nTot--;

    updateNeighbourLists(vacantNeighbours, neighbours, i, j, k);

}

void Solver::updateNeighbourLists(field<field<umat>> & A, field<field<umat>> & B, uint i, uint j, uint k)
{
    urowvec newRow = {i, j, k};

    //All comments works mirrored if we remove a particle

    //loop over the box surrounding the new site and
    //add the site as neighbours to surrounding sites
    //and then also remove it from the vacancy lists
    for (uint ci = 0; ci < 3; ++ci) {
        uint inp = (i + delta(ci) + N)%N;

        for (uint cj = 0; cj < 3; ++cj) {
            uint jnp = (j + delta(cj) + M)%M;

            for (uint ck = 0; ck < 3; ++ck) {

                if((ci == 1) && (cj == 1) && (ck == 1)) {
                    continue;
                }

                uint knp = (k + delta(ck) + L)%L;

                //                if ((abs((int)(k - knp))%L != 1) && (abs((int)(k - knp))%L != 29) &&(abs((int)(j - jnp))%M != 1)&& (abs((int)(j - jnp))%M != 29)  &&(abs((int)(i - inp))%N != 1) && (abs((int)(i - inp))%N != 29) ){
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
}
