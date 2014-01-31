#include "testbed.h"

#include <unittest++/UnitTest++.h>
#include <kMC>

#include <libconfig.h++>
#include <iostream>

using namespace libconfig;

testBed::testBed()
{
    Config cfg;

    cfg.readFile("infiles/config.cfg");

    const Setting & root = cfg.getRoot();

    solver = new KMCSolver(root);

    NX = solver->NX;
    NY = solver->NY;
    NZ = solver->NZ;


}

void testBed::testNeighbors()
{

    using namespace std;

    Site* site;
    int a;
    reset();

    for (uint x = 0; x < NX; ++x) {
        for (uint y = 0; y < NY; ++y) {
            for (uint z = 0; z < NZ; ++z) {

                site = solver->sites[x][y][z];
                for(uint i = 0; i < Site::neighborhoodLength(); ++i) {
                    for (uint j = 0; j < Site::neighborhoodLength(); ++j) {
                        for (uint k = 0; k < Site::neighborhoodLength(); ++k) {


                            if (i == Site::nNeighborsLimit() && j == Site::nNeighborsLimit() && k == Site::nNeighborsLimit()) {
                                CHECK_EQUAL(site->neighborHood[i][j][k], site);
                                continue;
                            }

                            nTrials++;
                            int delta = (int)(site->x()) - (int)(site->neighborHood[i][j][k]->x());

                            if ((delta + Site::originTransformVector()(i) + NX)%NX != 0) {

                                cout << "fail detected:" << endl;
                                cout << i << "  "<< x << "  " << site->neighborHood[i][j][k]->x() << endl;
                                cout << delta << "  " << Site::originTransformVector()(i) << "  " << NX << endl;
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


    CHECK_CLOSE(0.5, U, 0.01);
    CHECK_CLOSE(1.0/sqrt(12), stdU, 0.01);

    CHECK_CLOSE(0, N, 0.01);
    CHECK_CLOSE(1, stdN, 0.01);
}

void testBed::testBinarySearchChoise(uint LIM)
{
    uint choice;
    uint secondChoice;
    double R;

    solver->initialize();

    reset();

    solver->getRateVariables();

    while(nTrials < LIM)
    {

        R = solver->kTot*KMC_RNG_UNIFORM();

        choice = solver->getReactionChoice(R);
        secondChoice = 0;
        while (solver->accuAllRates.at(secondChoice) <= R) {
            secondChoice++;
        }

        secondChoice--;

        if (secondChoice == choice) {
            winCount++;
        } else {
            cout << choice << "  " << secondChoice << "  " << (int)choice - (int)secondChoice << endl;
            cout << solver->accuAllRates.at(choice) << "  " << solver->accuAllRates.at(secondChoice) << "  " << R << endl;
            cout << "--------" << endl;
            failCount++;
        }

        nTrials++;

    }

    CHECK_EQUAL(nTrials, winCount);

}

void testBed::testReactionChoise(uint LIM)
{
    uint choice;
    int old_count = -1;
    double kTot = 0;
    Reaction* reaction;

    solver->initialize();

    reset();

    solver->getRateVariables();

    while(nTrials < LIM)
    {

        for (double R : solver->accuAllRates) {

            choice = solver->getReactionChoice(R);

            kTot = 0;
            uint count = 0;
            for (uint i = 0; i < NX; ++i) {
                for (uint j = 0; j < NY; ++j) {
                    for (uint k = 0; k < NZ; ++k) {
                        for (Reaction* r : solver->sites[i][j][k]->activeReactions()) {

                            assert(r == solver->allReactions.at(count));

                            kTot += r->rate();

                            assert(kTot == solver->accuAllRates.at(count));

                            //if by adding this reaction we surpass the limit, we
                            //are done searching.
                            if (kTot >= R) {
                                reaction = r;
                                i = NX;
                                j = NY;
                                k = NZ;
                                break;
                            }

                            count++;
                        }
                    }
                }
            }

            if (!(choice == count)) {
                cout << "fail" << endl;
                cout << choice << " " << count << " " << (int)choice - (int)count << endl;
            }
            if (!(count == old_count + 1)) {
                cout << count << " " << solver->accuAllRates.size() << endl;
            }

            if (solver->allReactions[choice] == reaction) {
                winCount++;
            } else {
                failCount++;
            }

            nTrials++;
            old_count = count;
        }

    }

    CHECK_EQUAL(nTrials, winCount);

}
