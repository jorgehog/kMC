#include "testbed.h"

#include <unittest++/UnitTest++.h>
#include <kMC>

testBed::testBed()
{
}

void testBed::testNeighboursUpdating()
{

    uint NX = 10;
    uint NY = 10;
    uint NZ = 10;
    KMCSolver solver(NX, NY, NZ);


    solver.initialize();

    solver.getAllNeighbours();

    solver.getRateVariables();

    solver.allReactions[(uint)(KMC_RNG_UNIFORM()*solver.allReactions.size())]->execute();

    field<field<umat>> optneigh = solver.neighbours;

    solver.getAllNeighbours();

    uint failSize = 0;
    uint failVal = 0;
    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {

                //Check if we have the same number of neighbors, if not, failSize is incremented.
                if (solver.neighbours(i, j)(k).n_rows != optneigh(i, j)(k).n_rows) {
                    failSize++;
                }

                for (uint l = 0; l < solver.neighbours(i, j)(k).n_rows; ++l) {
                    for (int m = 0; m < 3; ++m) {

                        if (solver.neighbours(i, j)(k)(l, m) != optneigh(i, j)(k)(l, m)) {asd
                            failVal++;
                        }

                    }
                }
            }
        }
    }

    CHECK_EQUAL(0, n);

}

void testBed::testNeighbours()
{
    uint NX = 5;
    uint NY = 5;
    uint NZ = 5;
    KMCSolver solver(NX, NY, NZ);

    solver.nTot = 0;

    solver.sites[2][2][2]->activate();

    solver.sites[3][3][3]->activate();
    solver.sites[1][1][1]->activate();
    solver.sites[1][3][1]->activate();

    solver.getAllNeighbours();

    uint n = 0;



    solver.nTot++;


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
