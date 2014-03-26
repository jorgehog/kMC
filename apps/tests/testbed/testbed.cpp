#include "testbed.h"

#include "../snapshot/snapshot.h"

#include <unittest++/UnitTest++.h>

#include <iostream>

void testBed::makeSolver()
{

    solver = new KMCSolver();

    solver->setRNGSeed();

    DiffusionReaction::setPotentialParameters(1.0, 0.5, false);

    Site::setInitialBoundaries(Boundary::Periodic);

    Site::setInitialNNeighborsLimit(2, false);

    DiffusionReaction::setSeparation(1);

    Reaction::setBeta(0.5);

    solver->setBoxSize({10, 10, 10});


}

void testBed::testTotalParticleStateCounters()
{

    Site::resetBoundariesTo(Boundary::Periodic);

    CHECK_EQUAL(0, Site::totalActiveParticles(ParticleStates::surface));
    CHECK_EQUAL(0, accu(Site::totalActiveParticlesVector()));
    CHECK_EQUAL(NX()*NY()*NZ(), accu(Site::totalDeactiveParticlesVector()));
    CHECK_EQUAL(NX()*NY()*NZ(), Site::totalDeactiveParticles(ParticleStates::solution));

    initSimpleSystemParameters();

    activateAllSites();

    CHECK_EQUAL(0, Site::totalActiveParticles(ParticleStates::surface));
    CHECK_EQUAL(0, accu(Site::totalDeactiveParticlesVector()));
    CHECK_EQUAL(NX()*NY()*NZ(), accu(Site::totalActiveParticlesVector()));
    CHECK_EQUAL(NX()*NY()*NZ(), Site::totalActiveParticles(ParticleStates::solution));

    Site * boxCenter = solver->getSite(NX()/2, NY()/2, NZ()/2);

    boxCenter->deactivate();

    boxCenter->spawnAsFixedCrystal();

    CHECK_EQUAL(0, Site::totalActiveParticles(ParticleStates::surface));
    CHECK_EQUAL(0, accu(Site::totalDeactiveParticlesVector()));
    CHECK_EQUAL(NX()*NY()*NZ(), accu(Site::totalActiveParticlesVector()));
    CHECK_EQUAL(NX()*NY()*NZ() - 1, Site::totalActiveParticles(ParticleStates::crystal));

    boxCenter->deactivate();

    boxCenter->activate();

    CHECK_EQUAL(0, Site::totalActiveParticles(ParticleStates::surface));
    CHECK_EQUAL(0, accu(Site::totalDeactiveParticlesVector()));
    CHECK_EQUAL(NX()*NY()*NZ(), accu(Site::totalActiveParticlesVector()));
    CHECK_EQUAL(NX()*NY()*NZ(), Site::totalActiveParticles(ParticleStates::crystal));

    deactivateAllSites();

    CHECK_EQUAL(0, Site::totalActiveParticles(ParticleStates::surface));
    CHECK_EQUAL(0, accu(Site::totalActiveParticlesVector()));
    CHECK_EQUAL(NX()*NY()*NZ(), accu(Site::totalDeactiveParticlesVector()));
    CHECK_EQUAL(NX()*NY()*NZ(), Site::totalDeactiveParticles(ParticleStates::solution));

    boxCenter->spawnAsFixedCrystal();

    CHECK_EQUAL(0, Site::totalActiveParticles(ParticleStates::surface));
    CHECK_EQUAL(1, accu(Site::totalActiveParticlesVector()));
    CHECK_EQUAL(1, Site::totalActiveParticles(ParticleStates::fixedCrystal));
    CHECK_EQUAL(pow(DiffusionReaction::separation()*2 + 1, 3) - 1, Site::totalDeactiveParticles(ParticleStates::surface));

    solver->reset();

    CHECK_EQUAL(0, Site::totalActiveParticles(ParticleStates::surface));
    CHECK_EQUAL(0, accu(Site::totalActiveParticlesVector()));
    CHECK_EQUAL(NX()*NY()*NZ(), accu(Site::totalDeactiveParticlesVector()));
    CHECK_EQUAL(NX()*NY()*NZ(), Site::totalDeactiveParticles(ParticleStates::solution));


    uint C0 = NX()*NY()*NZ();
    uint C1 = 0;

    solver->forEachSiteDo([&C0, &C1] (Site * currentSite)
    {
        CHECK_EQUAL(C0, accu(Site::totalDeactiveParticlesVector()));
        CHECK_EQUAL(C0, Site::totalDeactiveParticles(ParticleStates::solution));

        CHECK_EQUAL(C1, accu(Site::totalActiveParticlesVector()));
        CHECK_EQUAL(C1, Site::totalActiveParticles(ParticleStates::solution));

        currentSite->activate();

        C0--;
        C1++;

    });

    deactivateAllSites();

    CHECK_EQUAL(0, Site::totalActiveParticles(ParticleStates::surface));
    CHECK_EQUAL(0, accu(Site::totalActiveParticlesVector()));
    CHECK_EQUAL(NX()*NY()*NZ(), accu(Site::totalDeactiveParticlesVector()));
    CHECK_EQUAL(NX()*NY()*NZ(), Site::totalDeactiveParticles(ParticleStates::solution));

    solver->setBoxSize({10, 10, 10});

    boxCenter = solver->getSite(NX()/2, NY()/2, NZ()/2);

    for (uint sep = 0; sep <= 3 ; ++sep)
    {

        Site::resetNNeighborsLimitTo(sep + 1, false);

        DiffusionReaction::resetSeparationTo(sep);

        boxCenter->spawnAsFixedCrystal();

        CHECK_EQUAL(pow(2*sep + 1, 3) - 1, Site::totalDeactiveParticles(ParticleStates::surface));

        boxCenter->deactivate();

    }

    DiffusionReaction::resetSeparationTo(3);

    solver->getSite(0, 0, 0)->spawnAsFixedCrystal();

    solver->forEachSiteDo_sendIndices([] (Site * site, uint x, uint y, uint z)
    {

        if ((x%2 == 0) && (y%2 == 0) && (z%2 == 0) && (!(x == 0 && y == 0 && z == 0)))
        {
            site->spawnAsFixedCrystal();
        }

    });

    solver->dumpXYZ();

    CHECK_EQUAL(NX()*NY()*NZ()/8, Site::totalActiveParticles(ParticleStates::fixedCrystal));
    CHECK_EQUAL(7*NX()*NY()*NZ()/8, Site::totalDeactiveParticles(ParticleStates::surface));



}

void testBed::testDistanceTo()
{

    int dx, dy, dz, dx2, dy2, dz2;
    uint adx, ady, adz;


    solver->setBoxSize({6, 6, 6}, false);
    Site::resetNNeighborsLimitTo(2);


    solver->forEachSiteDo([&] (Site * startSite)
    {
        solver->forEachSiteDo_sendIndices([&] (Site * endSite, uint endx, uint endy, uint endz)
        {

            Boundary::setupCurrentBoundaries(endx, endy, endz);

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

        });
    });

    initBoundaryTestParameters();

}

void testBed::testDeactivateSurface()
{

    solver->setBoxSize({10, 10, 10}, false);
    Site::resetNNeighborsToCrystallizeTo(1);

    Site * orig = solver->getSite(NX()/2, NY()/2, NZ()/2);
    Site * origNeighbor = solver->getSite(NX()/2+1, NY()/2, NZ()/2);
    Site * origNextNeighbor;
    Site * inBetweenSite;

    orig->spawnAsFixedCrystal();

    CHECK_EQUAL(ParticleStates::fixedCrystal, orig->particleState());
    CHECK_EQUAL(ParticleStates::surface, origNeighbor->particleState());


    uvec separations = {1, 2, 3};

    for (uint sep: separations)
    {

        Site::resetNNeighborsLimitTo(sep + 1);
        DiffusionReaction::resetSeparationTo(sep);

        orig->deactivate();
        orig->spawnAsFixedCrystal();


        CHECK_EQUAL(ParticleStates::fixedCrystal, orig->particleState());
        CHECK_EQUAL(ParticleStates::surface, origNeighbor->particleState());

        origNeighbor->activate();

        CHECK_EQUAL(ParticleStates::crystal, origNeighbor->particleState());

        origNextNeighbor = solver->getSite(NX()/2 + 1 + DiffusionReaction::separation(), NY()/2, NZ()/2);

        CHECK_EQUAL(ParticleStates::surface, origNextNeighbor->particleState());

        origNextNeighbor->activate();

        for (uint i = 1; i < DiffusionReaction::separation(); ++i)
        {
            inBetweenSite = solver->getSite(NX()/2 + 1 + i, NY()/2, NZ()/2);
            inBetweenSite->activate();
        }

        CHECK_EQUAL(ParticleStates::crystal, origNextNeighbor->particleState());

        Site::clearAffectedSites();
        for (uint i = 1; i < DiffusionReaction::separation(); ++i)
        {
            inBetweenSite = solver->getSite(NX()/2 + 1 + i, NY()/2, NZ()/2);
            inBetweenSite->deactivate();
        }

        origNeighbor->deactivate();

        CHECK_EQUAL(ParticleStates::fixedCrystal, orig->particleState());
        CHECK_EQUAL(ParticleStates::surface, origNeighbor->particleState());

        CHECK_EQUAL(false, origNextNeighbor->hasNeighboring(ParticleStates::crystal, DiffusionReaction::separation()));
        CHECK_EQUAL(ParticleStates::solution, origNextNeighbor->particleState());

        origNextNeighbor->deactivate();

    }

}

void testBed::testDiffusionSiteMatrixSetup()
{

    DiffusionReaction * currentDiffReaction;

    int i, j, k;

    solver->forEachSiteDo_sendIndices([&] (Site * site, uint x, uint y, uint z)
    {

        const Site & currentSite = *site;


        Boundary::setupCurrentBoundaries(x, y, z);


        for (Reaction * r : currentSite.reactions())
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


    });

}

void testBed::testNeighbors()
{

    int dx, dy, dz;

    solver->forEachSiteDo([&] (Site * currentSite)
    {
        currentSite->forEachNeighborDo_sendIndices([&] (Site * neighbor, uint i, uint j, uint k)
        {
            currentSite->distanceTo(neighbor, dx, dy, dz);

            CHECK_EQUAL(Site::originTransformVector(i), dx);
            CHECK_EQUAL(Site::originTransformVector(j), dy);
            CHECK_EQUAL(Site::originTransformVector(k), dz);

        });
    });
}

void testBed::testPropertyCalculations()
{

    initSimpleSystemParameters();

    CHECK_EQUAL(0, Site::nSurfaces());
    CHECK_EQUAL(0, Site::nCrystals());
    CHECK_EQUAL(0, Site::nSolutionParticles());
    CHECK_EQUAL(0, Site::getCurrentSolutionDensity());


    Site * center = getBoxCenter();

    center->spawnAsFixedCrystal();

    CHECK_EQUAL(26, Site::nSurfaces());
    CHECK_EQUAL(1, Site::nCrystals());
    CHECK_EQUAL(0, Site::nSolutionParticles());

    CHECK_EQUAL(0, Site::getCurrentSolutionDensity());
    CHECK_EQUAL(1.0, NX()*NY()*NZ()*Site::getCurrentRelativeCrystalOccupancy());

    umat boxTop = Site::getCurrentCrystalBoxTopology();

    CHECK_EQUAL(NX()/2, boxTop(0, 0));
    CHECK_EQUAL(NY()/2, boxTop(1, 0));
    CHECK_EQUAL(NZ()/2, boxTop(2, 0));

    CHECK_EQUAL(NX()/2, boxTop(0, 1));
    CHECK_EQUAL(NY()/2, boxTop(1, 1));
    CHECK_EQUAL(NZ()/2, boxTop(2, 1));

    solver->getSite(0, 0, 0)->activate();
    solver->getSite(1, 0, 0)->activate();
    solver->getSite(0, 1, 0)->activate();


    CHECK_EQUAL(3, Site::nSolutionParticles());

    CHECK_EQUAL(3.0/(NX()*NY()*NZ()- 1), Site::getCurrentSolutionDensity());

    activateAllSites();


    CHECK_EQUAL(0, Site::nSurfaces());
    CHECK_EQUAL(NX()*NY()*NZ(), Site::nCrystals());
    CHECK_EQUAL(0, Site::nSolutionParticles());

    CHECK_EQUAL(0, Site::getCurrentSolutionDensity());
    CHECK_EQUAL(1.0, Site::getCurrentRelativeCrystalOccupancy());

    boxTop = Site::getCurrentCrystalBoxTopology();

    CHECK_EQUAL(0, boxTop(0, 0));
    CHECK_EQUAL(0, boxTop(1, 0));
    CHECK_EQUAL(0, boxTop(2, 0));

    CHECK_EQUAL(NX()-1, boxTop(0, 1));
    CHECK_EQUAL(NY()-1, boxTop(1, 1));
    CHECK_EQUAL(NZ()-1, boxTop(2, 1));


    CHECK_EQUAL(0, Site::getCurrentSolutionDensity());

    deactivateAllSites();

    center->spawnAsFixedCrystal();

    uint sx = NX()/2;
    uint sy = NY()/2;
    uint sz = NZ()/2;

    solver->getSite(sx, sy    , sz + 1)->activate();
    solver->getSite(sx, sy    , sz + 2)->activate();

    solver->getSite(sx, sy - 1, sz - 1)->activate();
    solver->getSite(sx, sy - 1, sz    )->activate();
    solver->getSite(sx, sy - 1, sz + 1)->activate();


    solver->getSite(sx, sy + 1, sz + 1)->activate();

    /*            x
     *          x x x
     *          x F
     *          x
     *
     */


    boxTop = Site::getCurrentCrystalBoxTopology();

    CHECK_EQUAL(sx,     boxTop(0, 0));
    CHECK_EQUAL(sy - 1, boxTop(1, 0));
    CHECK_EQUAL(sz - 1, boxTop(2, 0));

    CHECK_EQUAL(sx,     boxTop(0, 1));
    CHECK_EQUAL(sy + 1, boxTop(1, 1));
    CHECK_EQUAL(sz + 2, boxTop(2, 1));

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
    for (uint i = 0; i < COUNT; ++i)
    {
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

    solver->initializeCrystal(0.3);

    solver->getRateVariables();

    uint N = 10;
    uint n = 0;

    while(n != N)
    {

        R = solver->kTot()*KMC_RNG_UNIFORM();

        choice = solver->getReactionChoice(R);

        secondChoice = 0;
        while (solver->accuAllRates().at(secondChoice) <= R)
        {
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


    solver->initializeCrystal(0.3);

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
                        bool notSet = true;
                        solver->getSite(i, j, k)->forEachActiveReactionDo([&] (Reaction * r)
                        {
                            if (!notSet)
                            {
                                return;
                            }

                            CHECK_EQUAL(r, solver->allPossibleReactions().at(count));

                            kTot += r->rate();

                            CHECK_EQUAL(kTot, solver->accuAllRates().at(count));

                            //if by adding this reaction we surpass the limit, we
                            //are done searching.
                            if (kTot > R)
                            {
                                reaction = r;

                                i = NX();
                                j = NY();
                                k = NZ();

                                notSet = false;
                            }

                            if (notSet)
                            {
                                count++;
                            }
                        });
                    }
                }
            }

            CHECK_EQUAL(solver->allPossibleReactions().at(choice), reaction);
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

    solver->initializeCrystal(0.3);

    solver->getRateVariables();

    solver->forEachSiteDo([&] (Site * site)
    {
        site->forEachActiveReactionDo([&] (Reaction * r)
        {

            E = ((DiffusionReaction*)r)->lastUsedEnergy();

            Esp = ((DiffusionReaction*)r)->lastUsedEsp();

            r->forceUpdateFlag(Reaction::defaultUpdateFlag);

            r->calcRate();

            CHECK_EQUAL(E, ((DiffusionReaction*)r)->lastUsedEnergy());

            CHECK_EQUAL(Esp, ((DiffusionReaction*)r)->lastUsedEsp());

        });
    });
}

void testBed::testEnergyAndNeighborSetup()
{

    int dx, dy, dz;
    uint ldx, ldy, ldz;

    double E;
    uint C;

    solver->initializeCrystal(0.3);

    uvec nn(Site::nNeighborsLimit());


    solver->forEachSiteDo([&] (Site * currentSite)
    {

        E = 0;
        C = 0;
        nn.zeros();

        solver->forEachSiteDo([&] (Site * otherSite)
        {

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
                                nn(Site::getLevel(ldx, ldy, ldz))++;

                                E += DiffusionReaction::potential(Site::nNeighborsLimit() + ldx,
                                                                  Site::nNeighborsLimit() + ldy,
                                                                  Site::nNeighborsLimit() + ldz);
                            }

                            C++;
                        }

                        CHECK_EQUAL(otherSite, currentSite->neighborhood(Site::nNeighborsLimit() + dx,
                                                                         Site::nNeighborsLimit() + dy,
                                                                         Site::nNeighborsLimit() + dz));

                        CHECK_EQUAL(currentSite, otherSite->neighborhood(Site::nNeighborsLimit() - dx,
                                                                         Site::nNeighborsLimit() - dy,
                                                                         Site::nNeighborsLimit() - dz));
                    }
                }
            }
        });

        uint nNeighbors = 0;
        currentSite->forEachNeighborDo([&nNeighbors] (Site * neighbor)
        {
            (void) neighbor;
            nNeighbors++;
        });

        CHECK_EQUAL(nNeighbors, C);

        for (uint K = 0; K < Site::nNeighborsLimit(); ++K)
        {
            CHECK_EQUAL(nn(K), currentSite->nNeighbors(K));
        }

        CHECK_CLOSE(E, currentSite->energy(), 0.00001);

    });

}

void testBed::testUpdateNeigbors()
{

    activateAllSites();

    double eMax = accu(DiffusionReaction::potentialBox());

    double blockedE;
    uvec nBlocked(Site::nNeighborsLimit());

    double accuBlockedE = 0;

    solver->forEachSiteDo([&] (Site * currentSite)
    {
        blockedE = 0;
        nBlocked.zeros();

        for (uint ni = 0; ni < Site::neighborhoodLength(); ++ni)
        {
            for (uint nj = 0; nj < Site::neighborhoodLength(); ++nj)
            {
                for (uint nk = 0; nk < Site::neighborhoodLength(); ++nk)
                {
                    if (currentSite->neighborhood(ni, nj, nk) == NULL)
                    {
                        nBlocked(Site::levelMatrix(ni, nj, nk))++;
                        blockedE += DiffusionReaction::potential(ni, nj, nk);
                    }
                }
            }
        }


        for (uint K = 0; K < Site::nNeighborsLimit(); ++K)
        {
            CHECK_EQUAL(2*(12*(K+1)*(K+1) + 1), currentSite->nNeighbors(K) + nBlocked(K));
        }

        CHECK_CLOSE(eMax, currentSite->energy() + blockedE, 0.001);

        accuBlockedE += blockedE;

    });

    CHECK_EQUAL(NX()*NY()*NZ(), Site::totalActiveSites());
    CHECK_CLOSE(NX()*NY()*NZ()*eMax, Site::totalEnergy() + accuBlockedE, 0.001);

    deactivateAllSites();

    solver->forEachSiteDo([&] (Site * site)
    {

        for (uint K = 0; K < Site::nNeighborsLimit(); ++K)
        {
            CHECK_EQUAL(0, site->nNeighbors(K));
        }

        CHECK_CLOSE(0, site->energy(), 0.001);

    });

    CHECK_EQUAL(0, Site::totalActiveSites());
    CHECK_CLOSE(0, Site::totalEnergy(), 0.001);

}


void testBed::testHasCrystalNeighbor()
{

    Site::resetNNeighborsLimitTo(2, false);
    solver->setBoxSize({10, 10, 10});

    Site::resetBoundariesTo(Boundary::Edge);

    DiffusionReaction::resetSeparationTo(1);

    Site::resetNNeighborsToCrystallizeTo(1);

    //Spawn a seed in the middle of the box.
    solver->getSite(NX()/2, NY()/2, NZ()/2)->spawnAsFixedCrystal();
    Site* initCrystal = solver->getSite(NX()/2, NY()/2, NZ()/2);

    Site *neighbor;
    uint level;

    //First we build a shell around the seed a distance 3 away which is all filled with particles.
    for (int i = -3; i < 4; ++i)
    {

        for (int j = -3; j < 4; ++j)
        {

            for (int k = -3; k < 4; ++k)
            {

                if (Site::getLevel(abs(i), abs(j), abs(k)) == 2)
                {
                    solver->getSite(NX()/2 + i, NY()/2 + j, NZ()/2 + k)->activate();
                    CHECK_EQUAL(ParticleStates::solution, solver->getSite(NX()/2 + i, NY()/2 + j, NZ()/2 + k)->particleState());
                }

            }

        }

    }

    uint nReactions = 8; //eight corners are free to move.
    for (int i = -2; i < 3; ++i)
    {

        for (int j = -2; j < 3; ++j)
        {

            for (int k = -2; k < 3; ++k)
            {

                uint level = Site::getLevel(abs(i), abs(j), abs(k));
                if (level == 1)
                {
                    nReactions += solver->getSite(NX()/2 + i, NY()/2 + j, NZ()/2 + k)->nNeighbors();
                }

            }

        }

    }


    for (uint i = 0; i < Site::neighborhoodLength(); ++i)
    {

        for (uint j = 0; j < Site::neighborhoodLength(); ++j)
        {

            for (uint k = 0; k < Site::neighborhoodLength(); ++k)
            {


                neighbor = initCrystal->neighborhood(i, j, k);

                //Then we check weather the middle is actually a crystal
                if (neighbor == initCrystal)
                {
                    assert(i == j && j == k && k == Site::nNeighborsLimit());

                    CHECK_EQUAL(ParticleStates::fixedCrystal, neighbor->particleState());

                    //it should not have any crystal neighbors
                    CHECK_EQUAL(false, neighbor->hasNeighboring(ParticleStates::crystal, DiffusionReaction::separation()));
                    continue;
                }

                level = Site::levelMatrix(i, j, k);

                //The first layer should now be a surface, which should be unblocked with a crystal neighbor.
                if (level == 0)
                {
                    CHECK_EQUAL(ParticleStates::surface, neighbor->particleState());
                    CHECK_EQUAL(true, neighbor->hasNeighboring(ParticleStates::crystal, DiffusionReaction::separation()));
                }

                //The second layer should be blocked because of the shell at distance 3, should be standard solution particles
                //without a crystal neighbor.
                else if (level == 1)
                {
                    CHECK_EQUAL(ParticleStates::solution, neighbor->particleState());
                    CHECK_EQUAL(false, neighbor->hasNeighboring(ParticleStates::crystal, DiffusionReaction::separation()));
                }

            }
        }
    }

    //deactivating the seed should bring everything to solutions except init seed which is surface.
    initCrystal->deactivate();

    //we now activate all neighbors. This should not make anything crystals.
    for (uint i = 0; i < 3; ++i)
    {

        for (uint j = 0; j < 3; ++j)
        {

            for (uint k = 0; k < 3; ++k)
            {
                neighbor = initCrystal->neighborhood(Site::nNeighborsLimit() - 1 + i,
                                                     Site::nNeighborsLimit() - 1 + j,
                                                     Site::nNeighborsLimit() - 1 + k);

                if (neighbor != initCrystal)
                {
                    neighbor->activate();
                }

            }

        }

    }

    for (uint i = 0; i < 3; ++i)
    {

        for (uint j = 0; j < 3; ++j)
        {

            for (uint k = 0; k < 3; ++k)
            {
                neighbor = initCrystal->neighborhood(Site::nNeighborsLimit() - 1 + i,
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
    for (int i = -3; i < 4; ++i)
    {

        for (int j = -3; j < 4; ++j)
        {

            for (int k = -3; k < 4; ++k)
            {

                if (Site::getLevel(abs(i), abs(j), abs(k)) == 0)
                {
                    CHECK_EQUAL(ParticleStates::crystal, solver->getSite(NX()/2 + i, NY()/2 + j, NZ()/2 + k)->particleState());
                }
                else if (Site::getLevel(abs(i), abs(j), abs(k)) == 1)
                {
                    CHECK_EQUAL(ParticleStates::surface, solver->getSite(NX()/2 + i, NY()/2 + j, NZ()/2 + k)->particleState());
                }
                else if (Site::getLevel(abs(i), abs(j), abs(k)) == 2)
                {
                    CHECK_EQUAL(ParticleStates::solution, solver->getSite(NX()/2 + i, NY()/2 + j, NZ()/2 + k)->particleState());
                    solver->getSite(NX()/2 + i, NY()/2 + j, NZ()/2 + k)->forEachActiveReactionDo([&nActives] (Reaction * r)
                    {
                        (void) r;
                        nActives++;
                    });
                }

            }

        }

    }

    //The number of possible reactions on the level=2 rim should be 1310
    CHECK_EQUAL(nReactions, nActives);



}

void testBed::testInitializationOfCrystal()
{

    solver->initializeCrystal(0.3);

    solver->forEachSiteDo([&] (Site * currentSite)
    {
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

    });

}

void testBed::testInitialReactionSetup()
{

    KMCDebugger_Init();

    Site * neighbor;

    uint nBlocked;

    solver->forEachSiteDo([&] (Site * currentSite)
    {

        if (currentSite->isActive())
        {

            CHECK_EQUAL(ParticleStates::fixedCrystal, currentSite->particleState());

            CHECK_EQUAL(0, currentSite->reactions().size());

        }

        else
        {

            nBlocked = 0;


            for (uint x = 0; x < 3; ++x)
            {
                for (uint y = 0; y < 3; ++y)
                {
                    for (uint z = 0; z < 3; ++z)
                    {

                        neighbor = currentSite->neighborhood(Site::nNeighborsLimit() - 1 + x,
                                                             Site::nNeighborsLimit() - 1 + y,
                                                             Site::nNeighborsLimit() - 1 + z);

                        if (neighbor == NULL)
                        {
                            nBlocked++;
                        }

                    }
                }
            }

            CHECK_EQUAL(26, currentSite->reactions().size() + nBlocked);

        }


    });

    solver->initializeCrystal(0.3);
    Site::updateAffectedSites();

    std::vector<Reaction*> oldReactions;
    double totRate1 = 0;

    solver->forEachSiteDo([&] (Site * currentSite)
    {
        currentSite->forEachActiveReactionDo([&] (Reaction * r)
        {
            KMCDebugger_AssertBool(currentSite->isActive(),
                                   "DEACTIVE SITE SHOULD HAVE NO REACTIONS",
                                   currentSite->info());

            KMCDebugger_AssertBool(r->isAllowed(),
                                   "REACTION NOT DEACTIVATED PROPERLY:",
                                   r->info());

            oldReactions.push_back(r);
            totRate1 += r->rate();

        });
    });

    solver->forEachActiveSiteDo([] (Site * currentSite)
    {
        currentSite->forEachActiveReactionDo([] (Reaction * r)
        {
                r->forceUpdateFlag(Reaction::defaultUpdateFlag);
        });

        currentSite->updateReactions();
    });


    std::vector<Reaction*> reactions;
    double totRate2 = 0;

    solver->forEachSiteDo([&] (Site * currentSite)
    {
        currentSite->forEachActiveReactionDo([&] (Reaction * r)
        {
            reactions.push_back(r);
            totRate2 += r->rate();
        });
    });

    CHECK_CLOSE(totRate1, totRate2, 0.000000001);
    CHECK_EQUAL(oldReactions.size(), reactions.size());

    for (uint i = 0; i < oldReactions.size(); ++i)
    {
        CHECK_EQUAL(reactions.at(i), oldReactions.at(i));
    }



}


void testBed::testSequential()
{

    solver->reset();

    initBoundaryTestParameters();


    const SnapShot s0(solver);

    const SnapShot & s1 = *testSequentialCore();


    solver->reset();

    const SnapShot & s2 = *testSequentialCore();


    CHECK_EQUAL(s1, s2);



    delete solver;

    makeSolver();

    initBoundaryTestParameters();

    const SnapShot s00(solver);

    const SnapShot & s3 = *testSequentialCore();


    delete solver;

    makeSolver();

    initBoundaryTestParameters();


    const SnapShot & s4 = *testSequentialCore();


    CHECK_EQUAL(s3, s4);


    CHECK_EQUAL(s0, s00);


    CHECK_EQUAL(s2, s3);


}

const SnapShot * testBed::testSequentialCore()
{

    uint nc = 1000;

    solver->setNumberOfCycles(nc);
    solver->setCyclesPerOutput(nc + 1);
    solver->initializeCrystal(0.2);

    solver->mainloop();

    return new SnapShot(solver);

}

void testBed::initBoundaryTestParameters()
{

    solver->setBoxSize({10, 10, 10}, false);

    Site::resetBoundariesTo(lastBoundaries);

    Site::resetNNeighborsLimitTo(3);

    DiffusionReaction::resetSeparationTo(1);

    Site::resetNNeighborsToCrystallizeTo(1);

    solver->setRNGSeed(Seed::specific, baseSeed);

}

void testBed::initSimpleSystemParameters()
{

    solver->setBoxSize({6, 6, 6}, false);

    DiffusionReaction::resetSeparationTo(1);

    Site::resetNNeighborsLimitTo(2);

    Site::resetNNeighborsToCrystallizeTo(1);

    Site::resetBoundariesTo(Boundary::Periodic);


}

void testBed::activateAllSites()
{

    solver->forEachSiteDo([] (Site * currentSite)
    {
        if (!currentSite->isActive())
        {
            currentSite->activate();
        }
    });

}

void testBed::deactivateAllSites()
{

    solver->forEachSiteDo([] (Site * currentSite)
    {
        if (currentSite->isActive())
        {
            currentSite->deactivate();
        }
    });
}

Site *testBed::getBoxCenter(const uint dx, const uint dy, const uint dz)
{
    return solver->getSite(NX()/2 + dx, NY()/2 + dy, NZ()/2 + dz);
}

void testBed::testKnownCase()
{

    delete solver;

    Config cfg;

    cfg.readFile("infiles/knowncase.cfg");

    const Setting & root = cfg.getRoot();

    solver = new KMCSolver(root);

    Site::resetBoundariesTo(lastBoundaries);

    ConcentrationWall::setMaxEventsPrCycle(5);

    solver->initializeCrystal(getSetting<double>(root, {"Initialization", "RelativeSeedSize"}));

    bool make = false;

    ifstream o;

    ofstream o2;

    stringstream fullName;
    fullName << "knowncase" << lastBoundariesName << ".txt";

    if (make)
    {
        o2.open(fullName.str());
    }

    else
    {
        o.open("infiles/" + fullName.str());

        if (!o.good())
        {
            cout << "NO KNOWNCASE FILE EXIST:" << "infiles/" + fullName.str() << endl;
            return;
        }
    }

    solver->mainloop();

    string line;
    stringstream s;

    uint equal = 0;

    solver->forEachSiteDo([&] (Site * site)
    {
        if (make)
        {
            o2 << site->isActive() << endl;
        }

        else
        {

            getline(o,line);
            s << site->isActive();
            if(s.str().compare(line) == 0)
            {
                equal++;
            }
            s.str(string());

        }
    });

    if(make)
    {
        o2.close();
        cout << "FILE MADE SUCCESSFULLY" << endl;
    }

    else
    {
        o.close();
        CHECK_EQUAL(NX()*NY()*NZ(), equal);
    }

    //if we dont reset to default solver, next time sequential test is called,
    //it will have mismatch in parameters such as
    //temperature and saturation.
    delete solver;
    makeSolver();

    ConcentrationWall::setMaxEventsPrCycle(3);


}

void testBed::testBoxSizes()
{


    uvec N = {6, 10, 15};

    Site::resetNNeighborsLimitTo(2, false);

    uvec3 boxSize;
    set<Site*> allSites;


    uint nBlocked;

    for (uint i = 0; i < N.n_elem; ++i)
    {

        for (uint j = 0; j < N.n_elem; ++j)
        {

            for (uint k = 0; k < N.n_elem; ++k)
            {

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

                solver->forEachSiteDo([&] (Site * currentSite)
                {
                    nBlocked = 0;

                    for (uint X = 0; X < Site::neighborhoodLength(); ++X)
                    {
                        for (uint Y = 0; Y < Site::neighborhoodLength(); ++Y)
                        {
                            for (uint Z = 0; Z < Site::neighborhoodLength(); ++Z)
                            {
                                if (currentSite->neighborhood(X, Y, Z) == NULL)
                                {
                                    nBlocked++;
                                }
                            }
                        }
                    }


                    uint nNeighbors = 0;
                    currentSite->forEachNeighborDo([&nNeighbors] (Site * neighbor)
                    {
                        (void) neighbor;
                        nNeighbors++;
                    });

                    CHECK_EQUAL(pow(Site::neighborhoodLength(), 3) - 1, nNeighbors + nBlocked);

                    currentSite->forEachNeighborDo([&allSites] (Site * neighbor)
                    {
                        allSites.insert(neighbor);
                    });

                });

                CHECK_EQUAL(nx*ny*nz, allSites.size());

            }
        }
    }



}

void testBed::testnNeiborsLimit()
{


    set<Site*> allSites;

    uvec nNlims = {1, 2, 3};

    uvec3 boxSize = {10, 10, 10};


    Site::resetBoundariesTo(Boundary::Periodic);

    solver->setBoxSize(boxSize, false);


    for (uint nNlim : nNlims)
    {
        allSites.clear();
        Site::resetNNeighborsLimitTo(nNlim, false);

        solver->forEachSiteDo([&] (Site * currentSite)
        {
            currentSite->forEachNeighborDo([&allSites] (Site * neighbor)
            {
                allSites.insert(neighbor);
            });

            uint nNeighbors = 0;
            currentSite->forEachNeighborDo([&nNeighbors] (Site * neighbor)
            {
                (void) neighbor;
                nNeighbors++;
            });

            CHECK_EQUAL(pow(2*nNlim + 1, 3) - 1, nNeighbors);

        });

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

    DiffusionReaction::resetSeparationTo(1);

    initialSeedSite->spawnAsFixedCrystal();

    trialSite->activate();


    for (uint nnts : nntcs)
    {

        Site::resetNNeighborsToCrystallizeTo(nnts);
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

        Site::clearAffectedSites();

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


        Site::clearAffectedSites();

    }



}

void testBed::testDiffusionSeparation()
{

    solver->setBoxSize({15, 15, 15}, false);


    Site * neighbor;
    Site * destination;
    Site * origin = solver->getSite(NX()/2, NY()/2, NZ()/2);

    origin->activate();

    uvec separations = {0, 1, 2, 3, 4, 5};

    for (uint sep : separations)
    {

        Site::resetNNeighborsLimitTo(sep + 1, false);

        DiffusionReaction::resetSeparationTo(sep);

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

                CHECK_EQUAL(neighbor->nActiveReactions(), neighbor->reactions().size());

                neighbor->deactivate();

                destination = solver->getSite(NX()/2 + sep + 1, NY()/2, NZ()/2);
                destination->activate();

                if (sep == 0)
                {
                    CHECK_EQUAL(25, destination->nActiveReactions());
                }

                else
                {

                    CHECK_EQUAL(17, destination->nActiveReactions());


                    destination->deactivate();

                    destination = solver->getSite(NX()/2 + sep, NY()/2, NZ()/2);

                    destination->activate();

                    CHECK_EQUAL(9, destination->nActiveReactions());


                    if (sep > 1)
                    {

                        destination->deactivate();

                        destination = solver->getSite(NX()/2 + sep - 1, NY()/2, NZ()/2);

                        destination->activate();

                        CHECK_EQUAL(0, destination->nActiveReactions());

                    }

                }

                destination->deactivate();
            }


            Site::clearAffectedSites();

        }

    }

    DiffusionReaction::resetSeparationTo(1);

}

void testBed::testAllPossibleRatesStuff()
{

    Reaction::setLinearRateScale(1000);

    solver->initializeCrystal(0.2);

    double kTot;
    vector<double> accuAllRates;
    vector<Reaction*> allPossibleReactions;

    uint cycle = 1;
    uint nCycles = 500;

    while (cycle <= nCycles)
    {

        solver->getRateVariables();

        fill_rate_stuff(accuAllRates, allPossibleReactions, kTot);

        CHECK_CLOSE(kTot, solver->kTot(), 1E-5);


        double R = solver->kTot()*KMC_RNG_UNIFORM();

        uint choice = solver->getReactionChoice(R);

        Reaction * selectedReaction = solver->allPossibleReactions().at(choice);

        selectedReaction->execute();

        Site::updateBoundaries();

        cycle++;

    }


    Reaction::setLinearRateScale(1);


}

void testBed::testReactionVectorUpdate()
{
    double c;

    Site * center = getBoxCenter();

    Reaction::setLinearRateScale(1000);


    //Activating the center, this should add 26 reactions with unit rate.
    center->activate();

    Site::updateAffectedSites();

    CHECK_EQUAL(0,  solver->m_availableReactionSlots.size());
    CHECK_EQUAL(26, solver->m_allPossibleReactions2.size());
    CHECK_EQUAL(26, solver->m_accuAllRates2.size());

    c = 0;
    for (uint i = 0; i < solver->m_allPossibleReactions2.size(); ++i)
    {
        c += solver->m_allPossibleReactions2.at(i)->rate();
        CHECK_EQUAL(c, solver->m_accuAllRates2.at(i));
    }

    CHECK_EQUAL(solver->kTot(), *(solver->m_accuAllRates2.end()-1));

    //Deactivating it should set all 26 reactions as available to overwrite.
    center->deactivate();

    Site::updateAffectedSites();

    CHECK_EQUAL(26,  solver->m_availableReactionSlots.size());
    CHECK_EQUAL(26,  solver->m_allPossibleReactions2.size());
    CHECK_EQUAL(26,  solver->m_accuAllRates2.size());

    CHECK_EQUAL(0, solver->kTot());

    for (uint i = 0; i < solver->m_allPossibleReactions2.size(); ++i)
    {
        CHECK_EQUAL(0, solver->m_accuAllRates2.at(i));
    }

    //Activating again should reset back to original case. Not trivial because this
    //time the vacancy list is not empty.
    center->activate();

    Site::updateAffectedSites();

    CHECK_EQUAL(0,  solver->m_availableReactionSlots.size());
    CHECK_EQUAL(26, solver->m_allPossibleReactions2.size());
    CHECK_EQUAL(26, solver->m_accuAllRates2.size());

    CHECK_EQUAL(26*Reaction::linearRateScale(), solver->kTot());

    c = 0;
    for (uint i = 0; i < solver->m_allPossibleReactions2.size(); ++i)
    {
        c += solver->m_allPossibleReactions2.at(i)->rate();
        CHECK_EQUAL(c, solver->m_accuAllRates2.at(i));
    }


    //activating a fixed crystal next to center particle. Fixed crystal has no
    //reactions so it should only induce 1 blocked reaction.
    Site * neighbor = getBoxCenter(1);
    neighbor->spawnAsFixedCrystal();

    Site::updateAffectedSites();

    CHECK_EQUAL(1,  solver->m_availableReactionSlots.size());

    c = 0;
    for (uint i = 0; i < solver->m_allPossibleReactions2.size(); ++i)
    {
        if (i == solver->m_availableReactionSlots.at(0))
        {
            continue;
        }

        c += solver->m_allPossibleReactions2.at(i)->rate();
        CHECK_CLOSE(c, solver->m_accuAllRates2.at(i), 0.00001);
    }

    //activating a new particle independent of the others should now only induce 25 more spots,
    //since one is already vacant.

    Site * distantCousin = getBoxCenter(0, Site::nNeighborsLimit() + 1);

    CHECK_EQUAL(true, distantCousin->isLegalToSpawn());

    distantCousin->activate();

    Site::updateAffectedSites();

    CHECK_EQUAL(0,  solver->m_availableReactionSlots.size());
    CHECK_EQUAL(26+25, solver->m_allPossibleReactions2.size());
    CHECK_EQUAL(26+25, solver->m_accuAllRates2.size());

    c = 0;
    for (uint i = 0; i < solver->m_allPossibleReactions2.size(); ++i)
    {
        c += solver->m_allPossibleReactions2.at(i)->rate();
        CHECK_CLOSE(c, solver->m_accuAllRates2.at(i), 0.00001);
    }


    //deactivating the neighbor should create two decoupled systems
    neighbor->deactivate();

    Site::updateAffectedSites();

    CHECK_EQUAL(0,    solver->m_availableReactionSlots.size());
    CHECK_EQUAL(2*26, solver->m_allPossibleReactions2.size());
    CHECK_EQUAL(2*26, solver->m_accuAllRates2.size());

    CHECK_CLOSE(2*26*Reaction::linearRateScale(), solver->kTot(), 0.00001);

    c = 0;
    for (uint i = 0; i < solver->m_allPossibleReactions2.size(); ++i)
    {
        c += solver->m_allPossibleReactions2.at(i)->rate();
        CHECK_CLOSE(c, solver->m_accuAllRates2.at(i), 0.00001);
    }

    //removing both the particles should bring us back to initial state

    center->deactivate();
    distantCousin->deactivate();

    Site::updateAffectedSites();

    CHECK_EQUAL(2*26, solver->m_availableReactionSlots.size());
    CHECK_EQUAL(2*26, solver->m_allPossibleReactions2.size());
    CHECK_EQUAL(2*26, solver->m_accuAllRates2.size());

    CHECK_CLOSE(0, solver->kTot(), 0.00001);

    for (uint i = 0; i < solver->m_allPossibleReactions2.size(); ++i)
    {
        CHECK_CLOSE(0, solver->m_accuAllRates2.at(i), 0.00001);
    }


    Reaction::setLinearRateScale(1);

}

void testBed::initBoundarySuite(const umat & boundaries)
{

    uint sum = accu(boundaries);


    string name = "";
    switch (sum)
    {
    case 6*Boundary::Periodic:
        name = "Periodic";
        break;
    case 6*Boundary::Edge:
        name = "Edge";
        break;
    case 6*Boundary::Surface:
        name = "Surface";
        break;
    case 6*Boundary::ConcentrationWall:
        name = "ConcentrationWall";
        break;
    default:
        name = "Mixed";
        break;
    }

    baseSeed = Seed::initialSeed;
    lastBoundaries = boundaries;
    lastBoundariesName = name;

    initBoundaryTestParameters();

}

void testBed::mainloop_meat()
{

    solver->getRateVariables();

    double R = solver->kTot()*KMC_RNG_UNIFORM();

    uint choice = solver->getReactionChoice(R);

    Reaction * selectedReaction = solver->allPossibleReactions().at(choice);

    selectedReaction->execute();

    Site::updateBoundaries();

}

void testBed::fill_rate_stuff(vector<double> & accuAllRates, vector<Reaction*> & allPossibleReactions, double & kTot)
{

    kTot = 0;
    accuAllRates.clear();
    allPossibleReactions.clear();

    double minRate = std::numeric_limits<double>::max();
    solver->forEachSiteDo([&] (Site * site)
    {
        site->forEachActiveReactionDo([&] (Reaction * reaction)
        {

            KMCDebugger_Assert(reaction->rate(), !=, Reaction::UNSET_RATE, "Reaction rate should not be unset at this point.", reaction->getFinalizingDebugMessage());

            kTot += reaction->rate();

            if (reaction->rate() < minRate)
            {
                minRate = reaction->rate();
            }

            accuAllRates.push_back(kTot);

            allPossibleReactions.push_back(reaction);

        });
    });
}


KMCSolver * testBed::solver;

wall_clock testBed::timer;

seed_type testBed::baseSeed;


string testBed::lastBoundariesName;

umat testBed::lastBoundaries;
