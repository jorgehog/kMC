#include "testbed.h"

#include "snapshot.h"

#include <unittest++/UnitTest++.h>

#include <iostream>

void testBed::makeSolver()
{

    Config cfg;

    cfg.readFile("infiles/config.cfg");

    const Setting & root = cfg.getRoot();

    solver = new KMCSolver(root);

}

void testBed::testDistanceTo()
{

    int dx, dy, dz, dx2, dy2, dz2;
    uint adx, ady, adz;

    Site* startSite;
    Site* endSite;


    for (uint startx = 0; startx < NX(); ++startx)
    {
        for (uint starty = 0; starty < NY(); ++starty)
        {
            for (uint startz = 0; startz < NZ(); ++startz)
            {

                startSite = solver->getSite(startx, starty, startz);

                for (uint endx = 0; endx < NX(); ++endx)
                {
                    for (uint endy = 0; endy < NY(); ++endy)
                    {
                        for (uint endz = 0; endz < NZ(); ++endz)
                        {

                            Boundary::setupCurrentBoundaries(endx, endy, endz);


                            endSite = solver->getSite(endx, endy, endz);

                            endSite->distanceTo(startSite, dx, dy, dz, true);

                            adx = dx;
                            ady = dy;
                            adz = dz;

                            endSite->distanceTo(startSite, dx, dy, dz);

                            CHECK_EQUAL(adx, abs(dx));
                            CHECK_EQUAL(ady, abs(dy));
                            CHECK_EQUAL(adz, abs(dz));

                            startSite->distanceTo(endSite, dx2, dy2, dz2);

                            if (adx != NX()/2)
                            {
                                CHECK_EQUAL(dx, -dx2);
                            }
                            else
                            {
                                CHECK_EQUAL(adx, abs(dx2));
                            }

                            if (ady != NY()/2)
                            {
                                CHECK_EQUAL(dy, -dy2);
                            }
                            else
                            {
                                CHECK_EQUAL(ady, abs(dy2));
                            }

                            if (adz != NZ()/2)
                            {
                                CHECK_EQUAL(dz, -dz2);
                            }
                            else
                            {
                                CHECK_EQUAL(adz, abs(dz2));
                            }

                            CHECK_EQUAL(startSite->x(), Boundary::currentBoundaries(0)->transformCoordinate(endSite->x() + dx));
                            CHECK_EQUAL(startSite->y(), Boundary::currentBoundaries(1)->transformCoordinate(endSite->y() + dy));
                            CHECK_EQUAL(startSite->z(), Boundary::currentBoundaries(2)->transformCoordinate(endSite->z() + dz));

                            if (Site::boundaryTypes(0) == Boundary::Periodic)
                            {
                                int X = (int)startSite->x() - (int)endSite->x();

                                if (X < -(int)NX()/2)
                                {
                                    X += NX();
                                }
                                else if (X > (int)NX()/2)
                                {
                                    X -= NX();
                                }


                                if (abs(X) == (int)NX()/2)
                                {
                                    CHECK_EQUAL(abs(X), adx);
                                }
                                else
                                {
                                    CHECK_EQUAL(dx, X);
                                }

                            }

                            if (Site::boundaryTypes(1) == Boundary::Periodic)
                            {
                                int Y = (int)startSite->y() - (int)endSite->y();

                                if (Y < -(int)NY()/2)
                                {
                                    Y += NY();
                                }
                                else if (Y > (int)NY()/2)
                                {
                                    Y -= NY();
                                }


                                if (abs(Y) == (int)NY()/2)
                                {
                                    CHECK_EQUAL(abs(Y), ady);
                                }
                                else
                                {
                                    CHECK_EQUAL(dy, Y);
                                }

                            }

                            if (Site::boundaryTypes(2) == Boundary::Periodic)
                            {
                                int Z = (int)startSite->z() - (int)endSite->z();

                                if (Z < -(int)NZ()/2)
                                {
                                    Z += NZ();
                                }
                                else if (Z > (int)NZ()/2)
                                {
                                    Z -= NZ();
                                }


                                if (abs(Z) == (int)NZ()/2)
                                {
                                    CHECK_EQUAL(abs(Z), adz);
                                }
                                else
                                {
                                    CHECK_EQUAL(dz, Z);
                                }
                            }

                        }
                    }
                }

            }
        }
    }

}

void testBed::testDeactivateSurface()
{


    Site::setNNeighborsToCrystallize(1);

    Site * orig = solver->getSite(NX()/2, NY()/2, NZ()/2);
    Site * origNeighbor = solver->getSite(NX()/2+1, NY()/2, NZ()/2);
    Site * origNextNeighbor;
    Site * inBetweenSite;

    orig->spawnAsFixedCrystal();

    CHECK_EQUAL(ParticleStates::fixedCrystal, orig->particleState());
    CHECK_EQUAL(ParticleStates::surface, origNeighbor->particleState());


    uvec separations = {1, 2, 3};
    solver->dumpXYZ();

    for (uint sep: separations)
    {

        Site::setNNeighborsLimit(sep + 1);
        DiffusionReaction::setSeparation(sep);

        orig->stripFixedCrystalProperty();
        orig->deactivate();
        orig->spawnAsFixedCrystal();

        solver->dumpXYZ();

        CHECK_EQUAL(ParticleStates::fixedCrystal, orig->particleState());
        CHECK_EQUAL(ParticleStates::surface, origNeighbor->particleState());

        origNeighbor->activate();
        solver->dumpXYZ();

        CHECK_EQUAL(ParticleStates::crystal, origNeighbor->particleState());

        origNextNeighbor = solver->getSite(NX()/2 + 1 + DiffusionReaction::separation(), NY()/2, NZ()/2);

        CHECK_EQUAL(ParticleStates::surface, origNextNeighbor->particleState());

        origNextNeighbor->activate();
        solver->dumpXYZ();
        //FILL EVERYTHING BETWEEN
        for (uint i = 1; i < DiffusionReaction::separation(); ++i)
        {
            inBetweenSite = solver->getSite(NX()/2 + 1 + i, NY()/2, NZ()/2);
            inBetweenSite->activate();
            solver->dumpXYZ();
        }

        CHECK_EQUAL(ParticleStates::crystal, origNextNeighbor->particleState());

        //DEACTIVATE EVERYTHING BETWEEN
        for (uint i = 1; i < DiffusionReaction::separation(); ++i)
        {
            inBetweenSite = solver->getSite(NX()/2 + 1 + i, NY()/2, NZ()/2);
            inBetweenSite->deactivate();
            solver->dumpXYZ();
        }

        origNeighbor->deactivate();
        solver->dumpXYZ();

        CHECK_EQUAL(ParticleStates::fixedCrystal, orig->particleState());
        CHECK_EQUAL(ParticleStates::surface, origNeighbor->particleState());

        CHECK_EQUAL(false, origNextNeighbor->hasNeighboring(ParticleStates::crystal, DiffusionReaction::separation()));
        CHECK_EQUAL(ParticleStates::solution, origNextNeighbor->particleState());

        origNextNeighbor->deactivate();
        solver->dumpXYZ();

    }

    solver->dumpXYZ();


}

void testBed::testDiffusionSiteMatrixSetup()
{

    DiffusionReaction * currentDiffReaction;

    int i, j, k;

    for (uint x = 0; x < NX(); ++x)
    {
        for (uint y = 0; y < NY(); ++y)
        {
            for (uint z = 0; z < NZ(); ++z)
            {

                const Site & currentSite = *(solver->getSite(x, y, z));


                Boundary::setupCurrentBoundaries(x, y, z);


                for (Reaction * r : currentSite.siteReactions())
                {

                    currentDiffReaction = (DiffusionReaction*)r;

                    const Site & site = *(currentDiffReaction->getReactionSite());
                    const Site & dest  = *(currentDiffReaction->destinationSite());

                    CHECK_EQUAL(currentSite, site);

                    site.distanceTo(&dest, i, j, k);

                    uint xt = Boundary::currentBoundaries(0)->transformCoordinate(x + i);
                    uint yt = Boundary::currentBoundaries(1)->transformCoordinate(y + j);
                    uint zt = Boundary::currentBoundaries(2)->transformCoordinate(z + k);

                    const Site & dest2 = *(solver->getSite(xt, yt, zt));

                    CHECK_EQUAL(dest, dest2);
                }


            }
        }
    }

}

void testBed::testNeighbors()
{

    Site* currentSite, *neighbor;

    int dx, dy, dz;

    for (uint x = 0; x < NX(); ++x)
    {
        for (uint y = 0; y < NY(); ++y)
        {
            for (uint z = 0; z < NZ(); ++z)
            {

                currentSite = solver->getSite(x, y, z);


                for(uint i = 0; i < Site::neighborhoodLength(); ++i)
                {
                    for (uint j = 0; j < Site::neighborhoodLength(); ++j)
                    {
                        for (uint k = 0; k < Site::neighborhoodLength(); ++k)
                        {

                            neighbor = currentSite->neighborHood(i, j, k);

                            if (neighbor == NULL)
                            {
                                continue;
                            }

                            if ((i == Site::nNeighborsLimit())
                                    && (j == Site::nNeighborsLimit())
                                    && (k == Site::nNeighborsLimit()))
                            {
                                CHECK_EQUAL(neighbor, currentSite);
                                continue;
                            }

                            currentSite->distanceTo(neighbor, dx, dy, dz);

                            CHECK_EQUAL(Site::originTransformVector(i), dx);
                            CHECK_EQUAL(Site::originTransformVector(j), dy);
                            CHECK_EQUAL(Site::originTransformVector(k), dz);


                        }
                    }
                }
            }
        }
    }
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

    solver->setRNGSeed(Seed::specific, Seed::initialSeed);

    for (uint i = 0; i < 1000000; ++i)
    {
        setu.push_back(KMC_RNG_UNIFORM());
        setn.push_back(KMC_RNG_NORMAL());
    }


    solver->setRNGSeed(Seed::specific, Seed::initialSeed);

    for (uint i = 0; i < 1000000; ++i)
    {
        CHECK_EQUAL(KMC_RNG_UNIFORM(), setu.at(i));
        CHECK_EQUAL(KMC_RNG_NORMAL(), setn.at(i));
    }


}

void testBed::testBinarySearchChoise()
{


    uint choice;
    uint secondChoice;
    double R;

    solver->initializeCrystal();

    solver->getRateVariables();

    uint N = 10;
    uint n = 0;

    while(n != N)
    {

        R = solver->kTot()*KMC_RNG_UNIFORM();

        choice = solver->getReactionChoice(R);

        secondChoice = 0;
        while (solver->accuAllRates().at(secondChoice) <= R) {
            secondChoice++;
        }

        CHECK_EQUAL(secondChoice, choice);

        n++;

    }



}

void testBed::testReactionChoise()
{


    uint choice, count, count2;

    double kTot, r_pre;

    Reaction* reaction;


    uint N = 3;
    uint n = 0;


    solver->initializeCrystal();

    solver->getRateVariables();

    while(n != N)
    {

        r_pre = 0;

        count2 = 0;

        for (double r_i : solver->accuAllRates())
        {

            double R = (r_i + r_pre)/2;

            choice = solver->getReactionChoice(R);

            kTot = 0;

            count = 0;

            for (uint i = 0; i < NX(); ++i)
            {
                for (uint j = 0; j < NY(); ++j)
                {
                    for (uint k = 0; k < NZ(); ++k)
                    {
                        for (Reaction* r : solver->getSite(i, j, k)->activeReactions())
                        {

                            CHECK_EQUAL(r, solver->allReactions().at(count));

                            kTot += r->rate();

                            CHECK_EQUAL(kTot, solver->accuAllRates().at(count));

                            //if by adding this reaction we surpass the limit, we
                            //are done searching.
                            if (kTot > R) {
                                reaction = r;

                                i = NX();
                                j = NY();
                                k = NZ();

                                break;
                            }

                            count++;

                        }
                    }
                }
            }

            CHECK_EQUAL(solver->allReactions().at(choice), reaction);
            CHECK_EQUAL(choice, count);
            CHECK_EQUAL(count, count2);

            count2++;

            r_pre = r_i;
        }

        n++;

    }


}

void testBed::testRateCalculation()
{


    double E, Esp;

    solver->initializeCrystal();

    solver->getRateVariables();

    for (uint i = 0; i < NX(); ++i)
    {
        for (uint j = 0; j < NY(); ++j)
        {
            for (uint k = 0; k < NZ(); ++k)
            {

                for (Reaction* r : solver->getSite(i, j, k)->activeReactions()) {

                    E = ((DiffusionReaction*)r)->lastUsedEnergy();

                    Esp = ((DiffusionReaction*)r)->lastUsedEsp();

                    r->forceUpdateFlag(Reaction::defaultUpdateFlag);

                    r->calcRate();

                    CHECK_EQUAL(E, ((DiffusionReaction*)r)->lastUsedEnergy());

                    CHECK_EQUAL(Esp, ((DiffusionReaction*)r)->lastUsedEsp());

                }
            }
        }
    }


}

void testBed::testEnergyAndNeighborSetup()
{

    int dx, dy, dz;
    uint ldx, ldy, ldz;

    double E;
    uint C, nBlocked;

    const Site* currentSite, * otherSite;

    solver->initializeCrystal();
    solver->dumpXYZ();

    uvec nn(Site::nNeighborsLimit());


    for (uint i = 0; i < NX(); ++i)
    {
        for (uint j = 0; j < NY(); ++j)
        {
            for (uint k = 0; k < NZ(); ++k)
            {

                E = 0;
                C = 0;
                nn.zeros();

                currentSite = solver->getSite(i, j, k);

                for (uint is = 0; is < NX(); ++is)
                {
                    for (uint js = 0; js < NY(); ++js)
                    {
                        for (uint ks = 0; ks < NZ(); ++ks)
                        {

                            otherSite = solver->getSite(is, js, ks);

                            currentSite->distanceTo(otherSite, dx, dy, dz);

                            ldx = abs(dx);

                            if (ldx <= Site::nNeighborsLimit())
                            {

                                ldy = abs(dy);
                                if (ldy <= Site::nNeighborsLimit())
                                {

                                    ldz = abs(dz);
                                    if (ldz <= Site::nNeighborsLimit())
                                    {

                                        if (currentSite != otherSite)
                                        {
                                            if (otherSite->isActive())
                                            {
                                                nn(Site::findLevel(ldx, ldy, ldz))++;

                                                E += DiffusionReaction::potential(Site::nNeighborsLimit() + ldx,
                                                                                  Site::nNeighborsLimit() + ldy,
                                                                                  Site::nNeighborsLimit() + ldz);
                                            }

                                            C++;
                                        }

                                        CHECK_EQUAL(otherSite, currentSite->neighborHood(Site::nNeighborsLimit() + dx,
                                                                                         Site::nNeighborsLimit() + dy,
                                                                                         Site::nNeighborsLimit() + dz));

                                        CHECK_EQUAL(currentSite, otherSite->neighborHood(Site::nNeighborsLimit() - dx,
                                                                                         Site::nNeighborsLimit() - dy,
                                                                                         Site::nNeighborsLimit() - dz));
                                    }
                                }
                            }
                        }
                    }
                }

                nBlocked = 0;

                for (Site * site : currentSite->allNeighbors())
                {
                    if (site == NULL)
                    {
                        nBlocked++;
                    }
                }


                CHECK_EQUAL(pow(Site::neighborhoodLength(), 3) - 1, C + nBlocked);

                for (uint K = 0; K < Site::nNeighborsLimit(); ++K) {
                    CHECK_EQUAL(nn(K), currentSite->nNeighbors(K));
                }

                CHECK_CLOSE(E, currentSite->energy(), 0.00001);

            }
        }
    }

}

void testBed::testUpdateNeigbors()
{

    return;

    bool enabled = KMCDebugger_IsEnabled;
    KMCDebugger_SetEnabledTo(false);

    CHECK_EQUAL(0, Site::totalEnergy());
    CHECK_EQUAL(0, Site::totalActiveSites());

    for (uint i = 0; i < NX(); ++i)
    {
        for (uint j = 0; j < NY(); ++j)
        {
            for (uint k = 0; k < NZ(); ++k)
            {
                solver->getSite(i, j, k)->activate();
            }
        }
    }

    double eMax = accu(DiffusionReaction::potentialBox());

    for (uint i = 0; i < NX(); ++i)
    {
        for (uint j = 0; j < NY(); ++j)
        {
            for (uint k = 0; k < NZ(); ++k)
            {

                for (uint K = 0; K < Site::nNeighborsLimit(); ++K)
                {
                    CHECK_EQUAL(2*(12*(K+1)*(K+1) + 1), solver->getSite(i, j, k)->nNeighbors(K));
                }

                CHECK_CLOSE(eMax, solver->getSite(i, j, k)->energy(), 0.001);

            }
        }
    }

    CHECK_EQUAL(NX()*NY()*NZ(), Site::totalActiveSites());
    CHECK_CLOSE(NX()*NY()*NZ()*eMax, Site::totalEnergy(), 0.001);

    for (uint i = 0; i < NX(); ++i)
    {
        for (uint j = 0; j < NY(); ++j)
        {
            for (uint k = 0; k < NZ(); ++k)
            {
                solver->getSite(i, j, k)->deactivate();
            }
        }
    }

    for (uint i = 0; i < NX(); ++i) {
        for (uint j = 0; j < NY(); ++j) {
            for (uint k = 0; k < NZ(); ++k) {

                for (uint K = 0; K < Site::nNeighborsLimit(); ++K) {
                    CHECK_EQUAL(0, solver->getSite(i, j, k)->nNeighbors(K));
                }

                CHECK_CLOSE(0, solver->getSite(i, j, k)->energy(), 0.001);

            }
        }
    }

    CHECK_EQUAL(0, Site::totalActiveSites());
    CHECK_CLOSE(0, Site::totalEnergy(), 0.001);

    KMCDebugger_SetEnabledTo(enabled);

}


void testBed::testHasCrystalNeighbor()
{

    Site::setNNeighborsLimit(2);

    //Spawn a seed in the middle of the box.
    solver->getSite(NX()/2, NY()/2, NZ()/2)->spawnAsFixedCrystal();
    Site* initCrystal = solver->getSite(NX()/2, NY()/2, NZ()/2);

    Site *neighbor;
    uint level;

    //First we build a shell around the seed a distance 3 away which is all filled with particles.
    for (int i = -3; i < 4; ++i) {

        for (int j = -3; j < 4; ++j) {

            for (int k = -3; k < 4; ++k) {

                if (Site::findLevel(abs(i), abs(j), abs(k)) == 2)
                {
                    solver->getSite(NX()/2 + i, NY()/2 + j, NZ()/2 + k)->activate();
                    CHECK_EQUAL(ParticleStates::solution, solver->getSite(NX()/2 + i, NY()/2 + j, NZ()/2 + k)->particleState());
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
                    nReactions += solver->getSite(NX()/2 + i, NY()/2 + j, NZ()/2 + k)->nNeighbors();
                }

            }

        }

    }


    for (uint i = 0; i < Site::neighborhoodLength(); ++i) {

        for (uint j = 0; j < Site::neighborhoodLength(); ++j) {

            for (uint k = 0; k < Site::neighborhoodLength(); ++k) {


                neighbor = initCrystal->neighborHood(i, j, k);

                //Then we check weather the middle is actually a crystal
                if (neighbor == initCrystal) {
                    assert(i == j && j == k && k == Site::nNeighborsLimit());

                    CHECK_EQUAL(ParticleStates::fixedCrystal, neighbor->particleState());

                    //it should not have any crystal neighbors
                    CHECK_EQUAL(false, neighbor->hasNeighboring(ParticleStates::crystal, DiffusionReaction::separation()));
                    continue;
                }

                level = Site::levelMatrix(i, j, k);

                //The first layer should now be a surface, which should be unblocked with a crystal neighbor.
                if (level == 0) {
                    CHECK_EQUAL(ParticleStates::surface, neighbor->particleState());
                    CHECK_EQUAL(true, neighbor->hasNeighboring(ParticleStates::crystal, DiffusionReaction::separation()));
                }

                //The second layer should be blocked because of the shell at distance 3, should be standard solution particles
                //without a crystal neighbor.
                else if (level == 1) {
                    CHECK_EQUAL(ParticleStates::solution, neighbor->particleState());
                    CHECK_EQUAL(false, neighbor->hasNeighboring(ParticleStates::crystal, DiffusionReaction::separation()));
                }

            }
        }
    }

    //deactivating the seed should bring everything to solutions except init seed which is surface.
    initCrystal->stripFixedCrystalProperty();
    initCrystal->deactivate();

    //we now activate all neighbors. This should not make anything crystals.
    for (uint i = 0; i < 3; ++i) {

        for (uint j = 0; j < 3; ++j) {

            for (uint k = 0; k < 3; ++k) {
                neighbor = initCrystal->neighborHood(Site::nNeighborsLimit() - 1 + i, Site::nNeighborsLimit() - 1 + j, Site::nNeighborsLimit() - 1 + k);

                if (neighbor != initCrystal) {
                    neighbor->activate();
                }

            }

        }

    }

    for (uint i = 0; i < 3; ++i) {

        for (uint j = 0; j < 3; ++j) {

            for (uint k = 0; k < 3; ++k) {
                neighbor = initCrystal->neighborHood(Site::nNeighborsLimit() - 1 + i,
                                                     Site::nNeighborsLimit() - 1 + j,
                                                     Site::nNeighborsLimit() - 1 + k);


                CHECK_EQUAL(ParticleStates::solution, neighbor->particleState());

            }

        }

    }

    //activating the seed. Should make closest neighbors crystals.
    initCrystal->spawnAsFixedCrystal();
    Site::updateAffectedSites();

    uint nActives = 0;
    for (int i = -3; i < 4; ++i) {

        for (int j = -3; j < 4; ++j) {

            for (int k = -3; k < 4; ++k) {

                if (Site::findLevel(abs(i), abs(j), abs(k)) == 0)
                {
                    CHECK_EQUAL(ParticleStates::crystal, solver->getSite(NX()/2 + i, NY()/2 + j, NZ()/2 + k)->particleState());
                }
                else if (Site::findLevel(abs(i), abs(j), abs(k)) == 1)
                {
                    CHECK_EQUAL(ParticleStates::surface, solver->getSite(NX()/2 + i, NY()/2 + j, NZ()/2 + k)->particleState());
                }
                else if (Site::findLevel(abs(i), abs(j), abs(k)) == 2)
                {
                    CHECK_EQUAL(ParticleStates::solution, solver->getSite(NX()/2 + i, NY()/2 + j, NZ()/2 + k)->particleState());
                    nActives += solver->getSite(NX()/2 + i, NY()/2 + j, NZ()/2 + k)->activeReactions().size();
                }

            }

        }

    }

    //The number of possible reactions on the level=2 rim should be 1310
    CHECK_EQUAL(nReactions, nActives);



}

void testBed::testInitializationOfCrystal()
{

    solver->initializeCrystal();

    Site* currentSite;

    for (uint i = 0; i < NX(); ++i)
    {
        for (uint j = 0; j < NY(); ++j)
        {
            for (uint k = 0; k < NZ(); ++k)
            {

                currentSite = solver->getSite(i, j, k);

                switch (currentSite->particleState())
                {
                case ParticleStates::solution:

                    //After initialization, a solution particle should not be blocked in any direction.
                    if (currentSite->isActive())
                    {
                        CHECK_EQUAL(0, currentSite->nNeighbors());
                    }

                    break;

                case ParticleStates::surface:

                    CHECK_EQUAL(true, !currentSite->isActive());
                    CHECK_EQUAL(true, currentSite->hasNeighboring(ParticleStates::crystal, DiffusionReaction::separation()));
                    CHECK_EQUAL(true, currentSite->nNeighbors() > 0);

                    break;

                case ParticleStates::crystal:

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



    KMCDebugger_Init();

    for (uint i = 0; i < NX(); ++i) {
        for (uint j = 0; j < NY(); ++j) {
            for (uint k = 0; k < NZ(); ++k) {

                CHECK_EQUAL(solver->getSite(i, j, k)->siteReactions().size(), 26);

            }
        }
    }

    solver->initializeCrystal();
    Site::updateAffectedSites();

    std::vector<Reaction*> oldReactions;
    double totRate1 = 0;

    for (uint i = 0; i < NX(); ++i) {
        for (uint j = 0; j < NY(); ++j) {
            for (uint k = 0; k < NZ(); ++k) {

                for (Reaction* r : solver->getSite(i, j, k)->activeReactions())
                {
                    KMCDebugger_AssertBool(solver->getSite(i, j, k)->isActive(),
                                           "DEACTIVE SITE SHOULD HAVE NO REACTIONS",
                                           solver->getSite(i, j, k)->info());

                    KMCDebugger_AssertBool(r->isAllowed(),
                                           "REACTION NOT DEACTIVATED PROPERLY:",
                                           r->info());

                    oldReactions.push_back(r);
                    totRate1 += r->rate();
                }
            }
        }
    }



    Site* currentSite;
    for (uint i = 0; i < NX(); ++i) {
        for (uint j = 0; j < NY(); ++j) {
            for (uint k = 0; k < NZ(); ++k) {

                currentSite = solver->getSite(i, j, k);

                currentSite->updateReactions();
                currentSite->calculateRates();

            }
        }
    }


    std::vector<Reaction*> reactions;
    double totRate2 = 0;

    for (uint i = 0; i < NX(); ++i) {
        for (uint j = 0; j < NY(); ++j) {
            for (uint k = 0; k < NZ(); ++k) {

                for (Reaction* r : solver->getSite(i, j, k)->activeReactions())
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
            cout << oldReactions.at(i)->info() << endl;
            exit(1);
        }
    }



}

void testBed::testSequential()
{



    uint nc = 100;

    solver->setNumberOfCycles(nc);
    solver->run();
    SnapShot s1(solver);

    seed_type seed = Seed::initialSeed;

    delete solver;

    makeSolver();
    KMC_INIT_RNG(seed);

    solver->setNumberOfCycles(nc);
    solver->run();
    SnapShot s2(solver);


    CHECK_EQUAL(s1, s2);


}

void testBed::testKnownCase()
{

    Config cfg;

    cfg.readFile("infiles/knowncase.cfg");

    const Setting & root = cfg.getRoot();

    solver = new KMCSolver(root);

    bool make = false;

    ifstream o;

    ofstream o2;

    if (make)
    {
        o2.open("knowncase.txt");
    }

    else
    {
        o.open("infiles/knowncase.txt");

        if (!o.good())
        {
            cout << "NO KNOWNCASE FILE EXIST." << endl;
            return;
        }
    }

    solver->run();

    string line;
    stringstream s;

    uint equal = 0;

    for (uint i = 0; i < NX(); ++i) {
        for (uint j = 0; j < NY(); ++j) {
            for (uint k = 0; k < NZ(); ++k) {

                if (make)
                {
                    o2 << solver->getSite(i, j, k)->isActive() << endl;
                }

                else
                {

                    getline(o,line);
                    s << solver->getSite(i, j, k)->isActive();
                    if(s.str().compare(line) == 0)
                    {
                        equal++;
                    }
                    s.str(string());

                }
            }
        }
    }

    if(make)
    {
        o2.close();
        cout << "FILE MADE SUCCESSFULLY" << endl;
    }

    else
    {
        o.close();
        CHECK_EQUAL(NX()*NY()*NZ(), equal);
        if (equal != NX()*NY()*NZ())
        {
            KMCDebugger_DumpFullTrace();
            exit(1);
        }
    }



}

void testBed::testBoxSizes()
{


    uvec N = {6, 10, 15};

    Site::setNNeighborsLimit(2);

    uvec3 boxSize;
    set<Site*> allSites;

    Site* currentSite;

    for (uint i = 0; i < N.n_elem; ++i) {

        for (uint j = 0; j < N.n_elem; ++j) {

            for (uint k = 0; k < N.n_elem; ++k) {

                uint nx = N(i);
                uint ny = N(j);
                uint nz = N(k);

                allSites.clear();

                boxSize = {nx, ny, nz};

                solver->setBoxSize(boxSize);

                CHECK_EQUAL(nx, NX());
                CHECK_EQUAL(ny, NY());
                CHECK_EQUAL(nz, NZ());

                CHECK_EQUAL(nx, Boundary::NX());
                CHECK_EQUAL(ny, Boundary::NY());
                CHECK_EQUAL(nz, Boundary::NZ());

                CHECK_EQUAL(nx, Site::NX());
                CHECK_EQUAL(ny, Site::NY());
                CHECK_EQUAL(nz, Site::NZ());

                CHECK_EQUAL(nx, Reaction::NX());
                CHECK_EQUAL(ny, Reaction::NY());
                CHECK_EQUAL(nz, Reaction::NZ());

                for (uint x = 0; x < NX(); ++x)
                {
                    for (uint y = 0; y < NY(); ++y)
                    {
                        for (uint z = 0; z < NZ(); ++z)
                        {

                            currentSite = solver->getSite(x, y, z);

                            CHECK_EQUAL(pow(Site::neighborhoodLength(), 3) - 1, currentSite->allNeighbors().size());

                            for (Site * neighbor : currentSite->allNeighbors())
                            {
                                allSites.insert(neighbor);
                            }
                        }
                    }
                }

                CHECK_EQUAL(nx*ny*nz, allSites.size());

            }
        }
    }



}

void testBed::testnNeiborsLimit()
{


    set<Site*> allSites;

    Site* currentSite;

    uvec nNlims = {1, 2, 3};

    uvec3 boxSize = {10, 10, 10};


    Site::setBoundaries(zeros<umat>(3, 2) + Boundary::Periodic);

    solver->setBoxSize(boxSize);


    for (uint nNlim : nNlims)
    {
        allSites.clear();
        Site::setNNeighborsLimit(nNlim);

        for (uint x = 0; x < NX(); ++x)
        {
            for (uint y = 0; y < NY(); ++y)
            {
                for (uint z = 0; z < NZ(); ++z)
                {

                    currentSite = solver->getSite(x, y, z);

                    allSites.insert(currentSite->allNeighbors().begin(), currentSite->allNeighbors().end());

                    CHECK_EQUAL(pow(2*nNlim + 1, 3) - 1, currentSite->allNeighbors().size());

                }
            }
        }

        CHECK_EQUAL(NX()*NY()*NZ(), allSites.size());

    }



}

void testBed::testnNeighborsToCrystallize()
{


    uvec nntcs = {1, 2, 3, 4, 5, 6, 7};

    Site * crystallizingSite;

    Site * initialSeedSite = solver->getSite(NX()/2, NY()/2, NZ()/2);

    Site * trialSite = solver->getSite(NX()/2 + 2, NY()/2, NZ()/2);

    uint totCrystalNeighbors;

    DiffusionReaction::setSeparation(1);

    initialSeedSite->spawnAsFixedCrystal();

    trialSite->activate();


    for (uint nnts : nntcs)
    {

        Site::setNNeighborsToCrystallize(nnts);
        totCrystalNeighbors = 0;

        //Fill a 3x3 surface with crystals.
        for (int i = -1; i <= 1; ++i)
        {
            for (int j = -1; j <= 1; ++j)
            {

                crystallizingSite = solver->getSite(NX()/2 + 1, NY()/2 + i, NZ()/2 + j);

                crystallizingSite->activate();

                totCrystalNeighbors++;

                CHECK_EQUAL(ParticleStates::crystal, crystallizingSite->particleState());


                if (totCrystalNeighbors >= nnts)
                {
                    CHECK_EQUAL(ParticleStates::crystal, trialSite->particleState());
                }

                else
                {
                    CHECK_EQUAL(ParticleStates::solution, trialSite->particleState());
                }


            }
        }

        for (int i = -1; i <= 1; ++i)
        {
            for (int j = -1; j <= 1; ++j)
            {

                crystallizingSite = solver->getSite(NX()/2 + 1, NY()/2 + i, NZ()/2 + j);

                crystallizingSite->deactivate();

                totCrystalNeighbors--;

                CHECK_EQUAL(ParticleStates::surface, crystallizingSite->particleState());


                if (totCrystalNeighbors >= nnts)
                {
                    CHECK_EQUAL(ParticleStates::crystal, trialSite->particleState());
                }

                else
                {
                    CHECK_EQUAL(ParticleStates::solution, trialSite->particleState());
                }

            }
        }

    }



}

void testBed::testDiffusionSeparation()
{


    bool enabled = KMCDebugger_IsEnabled;
    KMCDebugger_SetEnabledTo(false);

    solver->setBoxSize({15, 15, 15});


    Site * neighbor;
    Site * destination;
    Site * origin = solver->getSite(NX()/2, NY()/2, NZ()/2);

    origin->activate();

    uvec separations = {0, 1, 2, 3, 4, 5};

    for (uint sep : separations)
    {
        Site::setNNeighborsLimit(sep + 1);
        DiffusionReaction::setSeparation(sep);

        for (uint i = 1; i <= sep + 2; ++i)
        {

            neighbor = solver->getSite(NX()/2 + i, NY()/2, NZ()/2);

            bool allowed = neighbor->isLegalToSpawn();

            //sites only allowed to spawn if no reactions are blocked.
            //reaction blocked for up to sep + 1
            CHECK_EQUAL(!allowed, i <= sep + 1);

            if (i == sep + 2)
            {
                neighbor->activate();

                neighbor->updateReactions();

                CHECK_EQUAL(neighbor->activeReactions().size(), neighbor->siteReactions().size());

                neighbor->deactivate();

                destination = solver->getSite(NX()/2 + sep + 1, NY()/2, NZ()/2);
                destination->activate();

                destination->updateReactions();

                if (sep == 0)
                {
                    CHECK_EQUAL(25, destination->activeReactions().size());
                }

                else
                {

                    CHECK_EQUAL(17, destination->activeReactions().size());


                    destination->deactivate();

                    destination = solver->getSite(NX()/2 + sep, NY()/2, NZ()/2);

                    destination->activate();

                    destination->updateReactions();

                    CHECK_EQUAL(9, destination->activeReactions().size());


                    if (sep > 1)
                    {

                        destination->deactivate();

                        destination = solver->getSite(NX()/2 + sep - 1, NY()/2, NZ()/2);

                        destination->activate();

                        destination->updateReactions();

                        CHECK_EQUAL(0, destination->activeReactions().size());

                    }

                }

                destination->deactivate();

            }

        }

    }

    KMCDebugger_SetEnabledTo(enabled);


}

void testBed::testRunAllBoundaryTests(const umat & boundaries)
{

    uint sum = accu(boundaries);


    string name = "";
    switch (sum) {
    case 6*Boundary::Periodic:
        name = "Periodic";
        break;
    case 6*Boundary::Edge:
        name = "Edge";
        break;
    case 6*Boundary::Surface:
        name = "Surface";
        break;
    default:
        name = "Mixed";
        break;
    }

    cout << ".. for boundarytype " << name << endl;

    Site::setBoundaries(boundaries);
    Site::setNNeighborsLimit(3);
    solver->setBoxSize({10, 10, 10});

    cout << "   Running test DistanceTo" << endl;
    testDistanceTo();

    solver->reset();
    cout << "   Running test InitializationOfCrystal" << endl;
    testInitializationOfCrystal();

    solver->reset();
    cout << "   Running test DiffusionMatrixSetup" << endl;
    testDiffusionSiteMatrixSetup();

    solver->reset();
    cout << "   Running test Neighbors" << endl;
    testNeighbors();

    solver->reset();
    cout << "   Running test UpdateNeighbors" << endl;
    testUpdateNeigbors();

    solver->reset();
    cout << "   Running test EnergyAndNeighborSetup" << endl;
    testEnergyAndNeighborSetup();



}


KMCSolver * testBed::solver;

wall_clock testBed::timer;


