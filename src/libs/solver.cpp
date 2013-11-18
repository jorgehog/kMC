#include "solver.h"
#include "RNG/kMCRNG.h"

#include <sys/time.h>
#include <armadillo>
using namespace arma;

typedef bool*** BoolCube;


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

    uint N = 30;
    uint M = 30;
    uint L = 30;

    BoolCube siteCube = makeBoolCube(N, M, L);

    int a = 0;
    for (uint i = 0; i < N; ++i) {
        for (uint j = 0; j < M; ++j) {
            for (uint k = 0; k < L; ++k) {
                if (KMC_RNG_UNIFORM() > 0.5) {
                    siteCube[i][j][k] = true;
                    a++;
                }
            }
        }
    }

    cout << a/(double) (N*M*L) << endl;

    double T = 1E3;
    double t = 0;
    double dt;

    uint i, j, k;
    uint ip, jp, kp;

    while(t < T) {

        i = (uint) N*KMC_RNG_UNIFORM();
        j = (uint) M*KMC_RNG_UNIFORM();
        k = (uint) L*KMC_RNG_UNIFORM();

        if (siteCube[i][j][k]) {
            ip = (i + (int)(-1 + 2*KMC_RNG_UNIFORM()) + N)%N;
            jp = (j + (int)(-1 + 2*KMC_RNG_UNIFORM()) + M)%M;
            kp = (k + (int)(-1 + 2*KMC_RNG_UNIFORM()) + L)%L;

            if (!siteCube[ip][jp][kp]) {
                siteCube[ip][jp][kp] = true;
                siteCube[i][j][k] = false;
            }

        }

        dt = 1;
        t += dt;

    }

}

void Solver::test()
{
    cout << "test" << endl;
}
