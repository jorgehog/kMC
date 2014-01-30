#include "testbed.h"

#include <unittest++/UnitTest++.h>
#include <kMC>

#include <iostream>

testBed::testBed()
{
}

void testBed::testNeighbors()
{
    using namespace std;

    uint NX = 20;
    uint NY = 20;
    uint NZ = 20;

    KMCSolver* solver = new KMCSolver(0, NX, NY, NZ);

    Site* site;
    int a;
    reset();

    for (uint x = 0; x < NX; ++x) {
        for (uint y = 0; y < NY; ++y) {
            for (uint z = 0; z < NZ; ++z) {

                site = solver->sites[x][y][z];
                for(uint i = 0; i < Site::neighborhoodLength; ++i) {
                    for (uint j = 0; j < Site::neighborhoodLength; ++j) {
                        for (uint k = 0; k < Site::neighborhoodLength; ++k) {


                            if (i == Site::nNeighborsLimit && j == Site::nNeighborsLimit && k == Site::nNeighborsLimit) {
                                CHECK_EQUAL(site->neighborHood[i][j][k], site);
                                continue;
                            }

                            nTrials++;
                            int delta = (int)(site->x()) - (int)(site->neighborHood[i][j][k]->x());

                            if ((delta + Site::originTransformVector(i) + NX)%NX != 0) {

                                cout << "fail detected:" << endl;
                                cout << i << "  "<< x << "  " << site->neighborHood[i][j][k]->x() << endl;
                                cout << delta << "  " << Site::originTransformVector(i) << "  " << NX << endl;
                                cout << "------" << endl;

                                 cin >> a;
                                failCount++;
                            } else {
                                winCount++;
                            }

                        }
                    }
                }
            }
        }
    }

    CHECK_EQUAL(nTrials, winCount);
    CHECK_EQUAL(0, failCount);

}

void testBed::testRNG()
{
    KMC_INIT_RNG(time(NULL));

    double U  = 0;
    double U2 = 0;
    double N  = 0;
    double N2 = 0;

    double Ui;
    double Ni;

    uint COUNT = 1000000;
    for (uint i = 0; i < COUNT; ++i) {
        Ui = KMC_RNG_UNIFORM();
        Ni = KMC_RNG_NORMAL();

        U  += Ui;
        U2 += Ui*Ui;

        N  += Ni;
        N2 += Ni*Ni;
    }

    double stdU;
    double stdN;

    U  /= COUNT;
    U2 /= COUNT;
    N  /= COUNT;
    N2 /= COUNT;


    stdU = sqrt(U2 - U*U);
    stdN = sqrt(N2 - N*N);

    //    cout << U << "  " << N << endl;
    //    cout << stdU << "  " << stdN << endl;

    CHECK_CLOSE(0.5, U, 0.01);
    CHECK_CLOSE(1.0/sqrt(12), stdU, 0.01);

    CHECK_CLOSE(0, N, 0.01);
    CHECK_CLOSE(1, stdN, 0.01);
}
