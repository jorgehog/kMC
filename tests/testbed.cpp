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

testBed::~testBed()
{
    delete solver;
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
    double kTot;
    double r_pre;
    Reaction* reaction;

    reset();
    solver->initialize();
    solver->getRateVariables();

    while(nTrials < LIM)
    {

        cout << nTrials << " / " << LIM << endl;

        r_pre = 0;
        uint count2 = 0;
        for (double r_i : solver->accuAllRates) {

            double R = (r_i + r_pre)/2;
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
                            if (kTot > R) {
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

            if (solver->allReactions.at(choice) == reaction && choice == count && count == count2) {
                winCount++;
            } else {
                cout << r_i << "  " << r_pre << endl;
                cout << "res: "
                     << count2 << "  "
                     << choice << "  "
                     << count  << "  "
                     << solver->allReactions.size() - 1 << endl;
                return;
                failCount++;
            }

            count2++;
            r_pre = r_i;
        }

        nTrials++;
        CHECK_EQUAL(0, failCount);
    }
    cout << failCount << endl;
    cout << LIM << endl;
}

void testBed::testRateCalculation () {

    reset();

    solver->initialize();
    solver->getRateVariables();

    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {
                for (Reaction* r : solver->sites[i][j][k]->activeReactions()) {

                    double RATE = r->rate();
                    r->calcRate();
                    CHECK_EQUAL(r->rate(), RATE);

                    if (r->rate() < 1E-6 || r->rate() > 1000) {
                        double ESP = ((DiffusionReaction*)r)->getSaddleEnergy();
                        double E = solver->sites[i][j][k]->getEnergy();
                        cout << "RATE MESSED UP: "
                             << E
                             << "  "
                             << ESP
                             << "  "
                             << r->rate()
                             << endl;
                        failCount++;
                    }

                }
            }
        }
    }

    CHECK_EQUAL(0, failCount);
}

void testBed::testEnergyAndNeighborSetup()
{
    reset();

    solver->initialize();

    const Site* thisSite;
    const Site* otherSite;
    uvec nn(Site::nNeighborsLimit());

    int dx, dy, dz;
    uint ldx, ldy, ldz;
    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {

                thisSite = solver->sites[i][j][k];
                double E = 0;
                nn.zeros();
                uint C = 0;

                for (uint is = 0; is < NX; ++is) {
                    for (uint js = 0; js < NY; ++js) {
                        for (uint ks = 0; ks < NZ; ++ks) {

                            otherSite = solver->sites[is][js][ks];

                            thisSite->distanceTo(otherSite, dx, dy, dz);
                            ldx = abs(dx);
                            //                            cout << i << " " << j << " " << k << endl << is << " " << js << " " << ks << endl << dx << " " << dy << " " << dz << endl << "------" << endl;

                            if (ldx <= Site::nNeighborsLimit()) {
                                ldy = abs(dy);
                                if (ldy <= Site::nNeighborsLimit()) {
                                    ldz = abs(dz);

                                    if (ldz <= Site::nNeighborsLimit()) {

                                        if (thisSite != otherSite) {
                                            if (otherSite->active()) {
                                                nn(Site::getLevel(ldx, ldy, ldz))++;
                                                E += DiffusionReaction::potential()(Site::nNeighborsLimit() + ldx, Site::nNeighborsLimit() + ldy, Site::nNeighborsLimit() + ldz);
                                            }
                                            C++;
                                        }

                                        if (otherSite != thisSite->getNeighborhood()[Site::nNeighborsLimit() + dx][Site::nNeighborsLimit() + dy][Site::nNeighborsLimit() + dz]) {
                                            //                                            cout << "fail neighbor" << endl;
                                            failCount++;
                                        }

                                        if (thisSite != otherSite->getNeighborhood()[Site::nNeighborsLimit() - dx][Site::nNeighborsLimit() - dy][Site::nNeighborsLimit() - dz]) {
                                            //                                            cout << "fail neighbor" << endl;
                                            failCount++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                CHECK_EQUAL(pow(Site::neighborhoodLength(), 3) - 1, C);


                for (uint K = 0; K < Site::nNeighborsLimit(); ++K) {
                    CHECK_EQUAL(nn(K), thisSite->nNeighbors(K));
                }

                CHECK_CLOSE(E, thisSite->getEnergy(), 0.00001);

            }

        }
    }
    CHECK_EQUAL(0, failCount);
}

void testBed::testUpdateNeigbors()
{
    CHECK_EQUAL(0, Site::totalEnergy());
    CHECK_EQUAL(0, Site::totalActiveSites());

    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {
                solver->sites[i][j][k]->activate();
            }
        }
    }

    double eMax = accu(DiffusionReaction::potential());

    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {

                for (uint K = 0; K < Site::nNeighborsLimit(); ++K) {
                    CHECK_EQUAL(2*(12*(K+1)*(K+1) + 1), solver->sites[i][j][k]->nNeighbors(K));
                }

                CHECK_CLOSE(eMax, solver->sites[i][j][k]->getEnergy(), 0.00001);

            }
        }
    }

    CHECK_EQUAL(NX*NY*NZ, Site::totalActiveSites());
    CHECK_CLOSE(NX*NY*NZ*eMax, Site::totalEnergy(), 0.00001);

    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {
                solver->sites[i][j][k]->deactivate();
            }
        }
    }

    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {

                for (uint K = 0; K < Site::nNeighborsLimit(); ++K) {
                    CHECK_EQUAL(0, solver->sites[i][j][k]->nNeighbors(K));
                }

                CHECK_CLOSE(0, solver->sites[i][j][k]->getEnergy(), 0.00001);

            }
        }
    }

    CHECK_EQUAL(0, Site::totalActiveSites());
    CHECK_CLOSE(0, Site::totalEnergy(), 0.00001);

}



