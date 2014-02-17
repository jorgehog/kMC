#include "testbed.h"

#include "snapshot.h"

#include <kMC>

#include <unittest++/UnitTest++.h>

#include <iostream>


testBed::testBed()
{

    solver = makeSolver();

    NX = solver->NX;
    NY = solver->NY;
    NZ = solver->NZ;


}

KMCSolver *testBed::makeSolver()
{
    Config cfg;

    cfg.readFile("infiles/config.cfg");

    const Setting & root = cfg.getRoot();

    return new KMCSolver(root);

}

testBed::~testBed()
{
    delete solver;
}

void testBed::testDistanceTo()
{
    int dx, dy, dz;
    uint adx, ady, adz;

    int dx2, dy2, dz2;

    Site* start;
    Site* end;

    ivec deltax(NX);
    ivec deltay(NY);
    ivec deltaz(NZ);
    for(uint i = 0; i < NX; ++i)
    {
        deltax(i) = i;
        if (i > NX/2)
        {
            deltax(i) = -(int)(NX - i);
        }
    }
    for(uint i = 0; i < NY; ++i)
    {
        deltay(i) = i;
        if (i > NY/2)
        {
            deltay(i) = -(int)(NY - i);
        }
    }
    for(uint i = 0; i < NZ; ++i)
    {
        deltaz(i) = i;
        if (i > NZ/2)
        {
            deltaz(i) = -(int)(NZ - i);
        }
    }

    for (uint startx = 0; startx < NX; ++startx) {
        for (uint starty = 0; starty < NY; ++starty) {
            for (uint startz = 0; startz < NZ; ++startz) {

                start = solver->sites[startx][starty][startz];

                for (uint endx = 0; endx < NX; ++endx) {
                    for (uint endy = 0; endy < NY; ++endy) {
                        for (uint endz = 0; endz < NZ; ++endz) {


                            end = solver->sites[endx][endy][endz];

                            end->distanceTo(start, dx, dy, dz, true);

                            adx = dx;
                            ady = dy;
                            adz = dz;

                            end->distanceTo(start, dx, dy, dz);

                            CHECK_EQUAL(adx, abs(dx));
                            CHECK_EQUAL(ady, abs(dy));
                            CHECK_EQUAL(adz, abs(dz));

                            start->distanceTo(end, dx2, dy2, dz2);

                            if ((uint)abs(dx) != NX/2)
                            {
                                CHECK_EQUAL(dx, -dx2);
                            }
                            else
                            {
                                CHECK_EQUAL(abs(dx), abs(dx2));
                            }
                            if ((uint)abs(dy) != NY/2)
                            {
                                CHECK_EQUAL(dy, -dy2);
                            }
                            else
                            {
                                CHECK_EQUAL(abs(dy), abs(dy2));
                            }
                            if ((uint)abs(dz) != NZ/2)
                            {
                                CHECK_EQUAL(dz, -dz2);
                            }
                            else
                            {
                                CHECK_EQUAL(abs(dz), abs(dz2));
                            }

                            CHECK_EQUAL(dx2, deltax(((int)endx - (int)startx + NX)%NX));

                            CHECK_EQUAL(dy2, deltay(((int)endy - (int)starty + NY)%NY));

                            CHECK_EQUAL(dz2, deltaz(((int)endz - (int)startz + NZ)%NZ));

                        }
                    }
                }

            }
        }
    }

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
                                CHECK_EQUAL(site->m_neighborHood[i][j][k], site);
                                continue;
                            }

                            nTrials++;
                            int deltax = (int)(site->x()) - (int)(site->m_neighborHood[i][j][k]->x());
                            int deltay = (int)(site->y()) - (int)(site->m_neighborHood[i][j][k]->y());
                            int deltaz = (int)(site->z()) - (int)(site->m_neighborHood[i][j][k]->z());

                            if ((deltax + Site::originTransformVector()(i) + NX)%NX != 0) {

                                cout << "fail detected x:" << endl;
                                cout << i << "  "<< x << "  " << site->m_neighborHood[i][j][k]->x() << endl;
                                cout << deltax << "  " << Site::originTransformVector()(i) << "  " << NX << endl;
                                cout << "------" << endl;

                                cin >> a;
                                failCount++;
                            }
                            else if ((deltay + Site::originTransformVector()(j) + NY)%NY != 0) {

                                cout << "fail detected y:" << endl;
                                cout << j << "  "<< y << "  " << site->m_neighborHood[i][j][k]->y() << endl;
                                cout << deltay << "  " << Site::originTransformVector()(j) << "  " << NY << endl;
                                cout << "------" << endl;

                                cin >> a;
                                failCount++;

                            }
                            else if ((deltaz + Site::originTransformVector()(k) + NZ)%NZ != 0) {

                                cout << "fail detected z:" << endl;
                                cout << k << "  "<< z << "  " << site->m_neighborHood[i][j][k]->z() << endl;
                                cout << deltaz << "  " << Site::originTransformVector()(k) << "  " << NZ << endl;
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

    vector<double> setn;
    vector<double> setu;

    KMC_INIT_RNG(Seed::initialSeed);

    for (uint i = 0; i < 1000000; ++i)
    {
        setu.push_back(KMC_RNG_UNIFORM());
        setn.push_back(KMC_RNG_NORMAL());
    }

    KMC_RESET_RNG();

    for (uint i = 0; i < 1000000; ++i)
    {
        CHECK_EQUAL(KMC_RNG_UNIFORM(), setu.at(i));
        CHECK_EQUAL(KMC_RNG_NORMAL(), setn.at(i));
    }

}

void testBed::testBinarySearchChoise(uint LIM)
{
    uint choice;
    uint secondChoice;
    double R;

    solver->initializeCrystal();

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
    Reaction* reaction = NULL;

    reset();
    solver->initializeCrystal();
    solver->getRateVariables();

    while(nTrials < LIM)
    {

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
                cout << setprecision(16) << scientific << r_i - r_pre << endl;
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
}

void testBed::testRateCalculation () {

    reset();

    solver->initializeCrystal();
    solver->getRateVariables();

    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {

                for (Reaction* r : solver->sites[i][j][k]->activeReactions()) {

                    //                    double RATE = r->rate();
                    double E = ((DiffusionReaction*)r)->lastUsedE;
                    double Esp = ((DiffusionReaction*)r)->lastUsedEsp;
                    r->calcRate();

                    CHECK_EQUAL(E, ((DiffusionReaction*)r)->lastUsedE);

                    CHECK_EQUAL(Esp, ((DiffusionReaction*)r)->lastUsedEsp);

                    //                    CHECK_EQUAL(r->rate(), RATE);

                }
            }
        }
    }
}

void testBed::testEnergyAndNeighborSetup()
{
    reset();

    solver->initializeCrystal();

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
                                            if (otherSite->isActive()) {
                                                nn(Site::findLevel(ldx, ldy, ldz))++;
                                                E += DiffusionReaction::potential()(Site::nNeighborsLimit() + ldx, Site::nNeighborsLimit() + ldy, Site::nNeighborsLimit() + ldz);
                                            }
                                            C++;
                                        }

                                        if (otherSite != thisSite->neighborHood()[Site::nNeighborsLimit() + dx][Site::nNeighborsLimit() + dy][Site::nNeighborsLimit() + dz]) {
                                            //                                            cout << "fail neighbor" << endl;
                                            failCount++;
                                        }

                                        if (thisSite != otherSite->neighborHood()[Site::nNeighborsLimit() - dx][Site::nNeighborsLimit() - dy][Site::nNeighborsLimit() - dz]) {
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

                CHECK_CLOSE(E, thisSite->energy(), 0.00001);

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

                CHECK_CLOSE(eMax, solver->sites[i][j][k]->energy(), 0.001);

            }
        }
    }

    CHECK_EQUAL(NX*NY*NZ, Site::totalActiveSites());
    CHECK_CLOSE(NX*NY*NZ*eMax, Site::totalEnergy(), 0.001);

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

                CHECK_CLOSE(0, solver->sites[i][j][k]->energy(), 0.001);

            }
        }
    }

    CHECK_EQUAL(0, Site::totalActiveSites());
    CHECK_CLOSE(0, Site::totalEnergy(), 0.001);

}


void testBed::testHasCrystalNeighbor()
{

    //Spawn a seed in the middle of the box.
    solver->sites[NX/2][NY/2][NZ/2]->spawnAsCrystal();
    Site* initCrystal = solver->sites[NX/2][NY/2][NZ/2];

    Site *neighbor;
    uint level;

    //First we build a shell around the seed a distance 3 away which is all filled with particles.
    for (int i = -3; i < 4; ++i) {

        for (int j = -3; j < 4; ++j) {

            for (int k = -3; k < 4; ++k) {

                if (Site::findLevel(abs(i), abs(j), abs(k)) == 2)
                {
                    solver->sites[NX/2 + i][NY/2 + j][NZ/2 + k]->activate();
                    CHECK_EQUAL(particleState::solution, solver->sites[NX/2 + i][NY/2 + j][NZ/2 + k]->particleState());
                }

            }

        }

    }

    uint nReactions = 8; //eight corners are free to move.
    for (int i = -2; i < 3; ++i) {

        for (int j = -2; j < 3; ++j) {

            for (int k = -2; k < 3; ++k) {

                uint level = Site::findLevel(abs(i), abs(j), abs(k));
                if (level == 1)
                {
                    nReactions += solver->sites[NX/2 + i][NY/2 + j][NZ/2 + k]->nNeighbors();
                }

            }

        }

    }


    for (uint i = 0; i < Site::m_neighborhoodLength; ++i) {

        for (uint j = 0; j < Site::m_neighborhoodLength; ++j) {

            for (uint k = 0; k < Site::m_neighborhoodLength; ++k) {


                neighbor = initCrystal->m_neighborHood[i][j][k];

                //Then we check weather the middle is actually a crystal
                if (neighbor == initCrystal) {
                    assert(i == j && j == k && k == Site::m_nNeighborsLimit);

                    CHECK_EQUAL(particleState::crystal, neighbor->particleState());

                    //it should not have any crystal neighbors
                    CHECK_EQUAL(false, neighbor->hasNeighboring(particleState::crystal));
                    continue;
                }

                level = Site::m_levelMatrix(i, j, k);

                //The first layer should now be a surface, which should be unblocked with a crystal neighbor.
                if (level == 0) {
                    CHECK_EQUAL(particleState::surface, neighbor->particleState());
                    CHECK_EQUAL(true, neighbor->hasNeighboring(particleState::crystal));
                }

                //The second layer should be blocked because of the shell at distance 3, should be standard solution particles
                //without a crystal neighbor.
                else if (level == 1) {
                    CHECK_EQUAL(particleState::solution, neighbor->particleState());
                    CHECK_EQUAL(false, neighbor->hasNeighboring(particleState::crystal));
                }

            }
        }
    }

    //deactivating the seed should bring everything to solutions except init seed which is surface.
    initCrystal->deactivate();

    //we now activate all neighbors. This should not make anything crystals.
    for (uint i = 0; i < 3; ++i) {

        for (uint j = 0; j < 3; ++j) {

            for (uint k = 0; k < 3; ++k) {
                neighbor = initCrystal->m_neighborHood[Site::nNeighborsLimit() - 1 + i][Site::nNeighborsLimit() - 1 + j][Site::nNeighborsLimit() - 1 + k];

                if (neighbor != initCrystal) {
                    neighbor->activate();
                }

            }

        }

    }

    for (uint i = 0; i < 3; ++i) {

        for (uint j = 0; j < 3; ++j) {

            for (uint k = 0; k < 3; ++k) {
                neighbor = initCrystal->m_neighborHood[Site::nNeighborsLimit() - 1 + i][Site::nNeighborsLimit() - 1 + j][Site::nNeighborsLimit() - 1 + k];

                if (neighbor == initCrystal) {
                    CHECK_EQUAL(particleState::surface, neighbor->particleState());
                }
                else
                {
                    CHECK_EQUAL(particleState::solution, neighbor->particleState());
                }
            }

        }

    }

    //activating the seed. Should make closest neighbors crystals.
    initCrystal->activate();

    uint nActives = 0;
    for (int i = -3; i < 4; ++i) {

        for (int j = -3; j < 4; ++j) {

            for (int k = -3; k < 4; ++k) {

                if (Site::findLevel(abs(i), abs(j), abs(k)) == 0)
                {
                    CHECK_EQUAL(particleState::crystal, solver->sites[NX/2 + i][NY/2 + j][NZ/2 + k]->particleState());
                }
                else if (Site::findLevel(abs(i), abs(j), abs(k)) == 1)
                {
                    CHECK_EQUAL(particleState::surface, solver->sites[NX/2 + i][NY/2 + j][NZ/2 + k]->particleState());
                }
                else if (Site::findLevel(abs(i), abs(j), abs(k)) == 2)
                {
                    CHECK_EQUAL(particleState::solution, solver->sites[NX/2 + i][NY/2 + j][NZ/2 + k]->particleState());
                    nActives += solver->sites[NX/2 + i][NY/2 + j][NZ/2 + k]->activeReactions().size();
                }

            }

        }

    }

    //The number of possible reactions on the level=2 rim should be 1310
    CHECK_EQUAL(nReactions, nActives);

    solver->dumpXYZ();
}

void testBed::testInitializationOfCrystal()
{
    solver->initializeCrystal();

    Site* currentSite;
    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {

                currentSite = solver->sites[i][j][k];

                switch (currentSite->particleState()) {
                case particleState::solution:

                    //After initialization, a solution particle should not be blocked in any direction.
                    if (currentSite->isActive())
                    {
                        CHECK_EQUAL(0, currentSite->nNeighbors());
                    }

                    break;

                case particleState::surface:

                    CHECK_EQUAL(true, !currentSite->isActive());
                    CHECK_EQUAL(true, currentSite->hasNeighboring(particleState::crystal));
                    CHECK_EQUAL(true, currentSite->nNeighbors() > 0);

                    break;

                case particleState::crystal:

                    CHECK_EQUAL(true, currentSite->isActive());

                default:
                    break;
                }

            }
        }
    }

}

void testBed::testInitialReactionSetup()
{

    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {

                CHECK_EQUAL(solver->sites[i][j][k]->siteReactions().size(), 26);

            }
        }
    }

    solver->initializeCrystal();

    std::vector<Reaction*> oldReactions;
    double totRate1 = 0;

    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {

                for (Reaction* r : solver->sites[i][j][k]->activeReactions())
                {
                    if (!solver->sites[i][j][k]->isActive())
                    {
                        cout << "DEACTIVE SITE SHOULD HAVE NO REACTIONS" << endl;
                        solver->sites[i][j][k]->dumpInfo();
                        exit(1);
                    }

                    if (!r->isNotBlocked())
                    {
                        cout << "REACTION NOT DEACTIVATED PROPERLY:" << endl;
                        r->dumpInfo();

                        exit(1);
                    }
                    oldReactions.push_back(r);
                    totRate1 += r->rate();
                }
            }
        }
    }



    Site* currentSite;
    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {

                currentSite = solver->sites[i][j][k];

                currentSite->updateReactions();
                currentSite->calculateRates();

            }
        }
    }


    std::vector<Reaction*> reactions;
    double totRate2 = 0;

    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {

                for (Reaction* r : solver->sites[i][j][k]->activeReactions())
                {
                    reactions.push_back(r);
                    totRate2 += r->rate();
                }
            }
        }
    }

    CHECK_CLOSE(totRate1, totRate2, 0.000000001);
    CHECK_EQUAL(oldReactions.size(), reactions.size());

    for (uint i = 0; i < oldReactions.size(); ++i)
    {
        uint id1 = reactions.at(i)->ID();
        uint id2 = oldReactions.at(i)->ID();
        CHECK_EQUAL(id1, id2);

        if (id1 != id2)
        {
            cout << "reaction " << id2 << " mismatched: " << endl;
            oldReactions.at(i)->dumpInfo();
            exit(1);
        }
    }


}

void testBed::testSequential()
{
    solver->nCycles = 10;
    solver->run();
    SnapShot s1(solver);

    delete solver;
    solver = makeSolver();

    solver->nCycles = 10;
    solver->run();
    SnapShot s2(solver);

    CHECK_EQUAL(s1, s2);

}

void testBed::testKnownCase()
{

    CHECK_EQUAL(30, NX);
    CHECK_EQUAL(30, NY);
    CHECK_EQUAL(30, NZ);
    CHECK_EQUAL(2, Site::nNeighborsLimit());
    CHECK_EQUAL(5.0, Reaction::beta);
    CHECK_EQUAL(1.0, Reaction::m_linearRateScale);
    CHECK_EQUAL(0.5, DiffusionReaction::rPower);
    CHECK_EQUAL(1.0, DiffusionReaction::scale);
    CHECK_EQUAL(0.3, solver->saturation);
    CHECK_EQUAL(0.3, solver->RelativeSeedSize);
    CHECK_EQUAL(10000, solver->nCycles);
    CHECK_EQUAL(1000, solver->cyclesPerOutput);
    CHECK_EQUAL(1392202630, Seed::initialSeed);

    solver->nCycles = 4;
    solver->run();

    ifstream o;
    o.open("stuff.txt");
    string line;
    stringstream s;

    reset();

    for (uint i = 0; i < NX; ++i) {
        for (uint j = 0; j < NY; ++j) {
            for (uint k = 0; k < NZ; ++k) {
                getline(o,line);
                s << solver->sites[i][j][k]->isActive();
                if(s.str().compare(line) == 0)
                {
                    winCount++;
                }
                s.str(string());
            }
        }
    }
    o.close();

    CHECK_EQUAL(NX*NY*NZ, winCount);


}



