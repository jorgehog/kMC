#include "testbed.h"

#include "../snapshot/snapshot.h"

#include "../dummyreaction.h"

#include <unittest++/UnitTest++.h>

#include <iostream>

void testBed::makeSolver()
{

    solver = new KMCSolver();

    solver->setRNGSeed();

    solver->setNumberOfCycles(1000);

    solver->setTargetConcentration(0.01);

    DiffusionReaction::setPotentialParameters(1.0, 0.5, false);

    Site::setInitialBoundaries(Boundary::Periodic);

    Site::setInitialNNeighborsLimit(2, false);

    Reaction::setBeta(0.5);

    solver->setBoxSize({10, 10, 10});


}

void testBed::testTotalParticleStateCounters()
{

    Site::resetBoundariesTo(Boundary::Periodic);

    CHECK_EQUAL(0, SoluteParticle::nParticles());
    CHECK_EQUAL(0, accu(SoluteParticle::totalParticlesVector()));

    initSimpleSystemParameters();

    activateAllSites();

    CHECK_EQUAL(0, SoluteParticle::totalParticles(ParticleStates::surface));
    CHECK_EQUAL(NX()*NY()*NZ(), accu(SoluteParticle::totalParticlesVector()));
    CHECK_EQUAL(NX()*NY()*NZ(), SoluteParticle::totalParticles(ParticleStates::crystal));

    deactivateAllSites();

    CHECK_EQUAL(0, SoluteParticle::totalParticles(ParticleStates::surface));
    CHECK_EQUAL(0, accu(SoluteParticle::totalParticlesVector()));

    //    boxCenter->spawnAsFixedCrystal();

    //    CHECK_EQUAL(0, Site::totalActiveParticles(ParticleStates::surface));
    //    CHECK_EQUAL(1, accu(Site::totalActiveParticlesVector()));
    //    CHECK_EQUAL(1, Site::totalActiveParticles(ParticleStates::fixedCrystal));
    //    CHECK_EQUAL(pow(DiffusionReaction::separation()*2 + 1, 3) - 1, Site::totalDeactiveParticles(ParticleStates::surface));

    solver->reset();

    CHECK_EQUAL(0, accu(SoluteParticle::totalParticlesVector()));

    uint C1 = 0;

    solver->forEachSiteDo([&C1] (Site * currentSite)
    {
        CHECK_EQUAL(C1, accu(SoluteParticle::totalParticlesVector()));

        solver->forceSpawnParticle(currentSite);
        C1++;
    });

    deactivateAllSites();

    CHECK_EQUAL(0, accu(SoluteParticle::totalParticlesVector()));

    solver->setBoxSize({10, 10, 10});

    solver->forEachSiteDo([] (Site * site)
    {
        if (!site->isActive())
        {
            return;
        }

        if ((site->x()%2 == 0) && (site->y()%2 == 0) && (site->z()%2 == 0))
        {
            solver->forceSpawnParticle(site);
        }

    });

    CHECK_EQUAL(NX()*NY()*NZ()/8, SoluteParticle::totalParticles(ParticleStates::solvant));



}

void testBed::testDistanceTo()
{

    int dx, dy, dz, dx2, dy2, dz2;
    uint adx, ady, adz;


    solver->setBoxSize({6, 6, 6}, false);
    Site::resetNNeighborsLimitTo(2);


    solver->forEachSiteDo([&] (Site * startSite)
    {
        solver->forEachSiteDo([&] (Site * endSite)
        {

            Boundary::setupCurrentBoundaries(endSite->x(), endSite->y(), endSite->z());

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


void testBed::testDiffusionSiteMatrixSetup()
{

    DiffusionReaction * currentDiffReaction;

    int _i, _j, _k;

    solver->initializeCrystal(0.2);

    for (SoluteParticle *particle : solver->particles())
    {

        Boundary::setupCurrentBoundaries(particle->x(), particle->y(), particle->z());

        for (uint i = 0; i < 3; ++i)
        {
            for (uint j = 0; j < 3; ++j)
            {
                for (uint k = 0; k < 3; ++k)
                {

                    if (i == j && j == k && k == 1)
                    {
                        continue;
                    }

                    currentDiffReaction = particle->diffusionReactions(i, j, k);

                    const Site & site = *(currentDiffReaction->site());
                    const Site & dest  = *(currentDiffReaction->destinationSite());

                    CHECK_EQUAL(particle->site(), &site);

                    site.distanceTo(&dest, _i, _j, _k);

                    uint xt = Boundary::currentBoundaries(0)->transformCoordinate(particle->x() + _i);
                    uint yt = Boundary::currentBoundaries(1)->transformCoordinate(particle->y() + _j);
                    uint zt = Boundary::currentBoundaries(2)->transformCoordinate(particle->z() + _k);

                    CHECK_EQUAL(i, _i+1);
                    CHECK_EQUAL(j, _j+1);
                    CHECK_EQUAL(k, _k+1);

                    const Site & dest2 = *(solver->getSite(xt, yt, zt));

                    CHECK_EQUAL(dest, dest2);


                }
            }
        }


    }

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

    CHECK_EQUAL(0, SoluteParticle::nSurfaces());
    CHECK_EQUAL(0, SoluteParticle::nCrystals());
    CHECK_EQUAL(0, SoluteParticle::nSolutionParticles());
    CHECK_EQUAL(0, SoluteParticle::getCurrentConcentration());

    solver->setTargetConcentration(0);

    double relativeSeedSize = 0.4;
    solver->initializeCrystal(relativeSeedSize);

    CHECK_EQUAL(0, SoluteParticle::getCurrentConcentration());
    CHECK_CLOSE(pow(relativeSeedSize, 3), SoluteParticle::getCurrentRelativeCrystalOccupancy(), 1E-5);
    solver->dumpXYZ(1337);

    umat boxTop = Site::getCurrentCrystalBoxTopology();

    uint crystalSizeX = round(NX()*relativeSeedSize);
    uint crystalSizeY = round(NY()*relativeSeedSize);
    uint crystalSizeZ = round(NZ()*relativeSeedSize);

    uint crystalStartX = NX()/2 - crystalSizeX/2;
    uint crystalStartY = NY()/2 - crystalSizeY/2;
    uint crystalStartZ = NZ()/2 - crystalSizeZ/2;

    uint crystalEndX = crystalStartX + crystalSizeX;
    uint crystalEndY = crystalStartY + crystalSizeY;
    uint crystalEndZ = crystalStartZ + crystalSizeZ;


    CHECK_EQUAL(crystalStartX+1, boxTop(0, 0));
    CHECK_EQUAL(crystalStartY+1, boxTop(1, 0));
    CHECK_EQUAL(crystalStartZ+1, boxTop(2, 0));

    CHECK_EQUAL(crystalEndX-2, boxTop(0, 1));
    CHECK_EQUAL(crystalEndY-2, boxTop(1, 1));
    CHECK_EQUAL(crystalEndZ-2, boxTop(2, 1));

    solver->forceSpawnParticle(solver->getSite(0, 0, 0));
    solver->forceSpawnParticle(solver->getSite(2, 0, 0));
    solver->forceSpawnParticle(solver->getSite(0, 2, 0));


    CHECK_EQUAL(3, SoluteParticle::nSolutionParticles());

    CHECK_EQUAL(3.0/(NX()*NY()*NZ()- crystalSizeX*crystalSizeY*crystalSizeZ), SoluteParticle::getCurrentConcentration());

    activateAllSites();


    CHECK_EQUAL(0, SoluteParticle::nSurfaces());
    CHECK_EQUAL(NX()*NY()*NZ(), SoluteParticle::nCrystals());
    CHECK_EQUAL(0, SoluteParticle::nSolutionParticles());

    CHECK_EQUAL(0, SoluteParticle::getCurrentConcentration());
    CHECK_EQUAL(1.0, SoluteParticle::getCurrentRelativeCrystalOccupancy());

    boxTop = Site::getCurrentCrystalBoxTopology();

    CHECK_EQUAL(0, boxTop(0, 0));
    CHECK_EQUAL(0, boxTop(1, 0));
    CHECK_EQUAL(0, boxTop(2, 0));

    CHECK_EQUAL(NX()-1, boxTop(0, 1));
    CHECK_EQUAL(NY()-1, boxTop(1, 1));
    CHECK_EQUAL(NZ()-1, boxTop(2, 1));



    deactivateAllSites();

    CHECK_EQUAL(0, SoluteParticle::getCurrentConcentration());

    CHECK_EQUAL(0, SoluteParticle::nParticles());


    uint sx = NX()/2;
    uint sy = NY()/2;
    uint sz = NZ()/2;

    vector<Site*> crystalSites = {solver->getSite(sx, sy    , sz    ),
                                  solver->getSite(sx, sy    , sz    ),

                                  solver->getSite(sx, sy    , sz + 1),
                                  solver->getSite(sx, sy    , sz + 2),

                                  solver->getSite(sx, sy - 1, sz - 1),
                                  solver->getSite(sx, sy - 1, sz    ),
                                  solver->getSite(sx, sy - 1, sz + 1),

                                  solver->getSite(sx, sy + 1, sz + 1)
                                 };

    uint x, y, z;

    Site *currentSite;

    for (Site * site : crystalSites)
    {
        x = site->x();
        y = site->y();
        z = site->z();

        for (int i = -1; i <= 1; ++i)
        {
            for (int j = -1; j <= 1; ++j)
            {
                for (int k = -1; k <= 1; ++k)
                {

                    currentSite = solver->getSite(i + x, j + y, k + z);

                    if (!currentSite->isActive())
                    {
                        solver->forceSpawnParticle(currentSite);
                    }
                }
            }
        }
    }

    /*            x
     *          x x x
     *          x x
     *          x
     *
     */

    solver->dumpXYZ(1338);

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

    KMCDebugger_Assert(0, ==, SoluteParticle::affectedParticles().size());
    KMCDebugger_Assert(0, ==, solver->allPossibleReactions().size());
    KMCDebugger_Assert(0, ==, solver->accuAllRates().size());

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

            for (SoluteParticle *particle : solver->particles())
            {

                bool notSet = true;

                particle->forEachActiveReactionDo([&] (Reaction *r)
                {

                    if (!notSet)
                    {
                        return;
                    }

                    CHECK_EQUAL(*r, *solver->allPossibleReactions().at(count));

                    kTot += r->rate();

                    CHECK_EQUAL(kTot, solver->accuAllRates().at(count));

                    //if by adding this reaction we surpass the limit, we
                    //are done searching.
                    if (kTot > R)
                    {
                        reaction = r;

                        notSet = false;

                        return;
                    }

                    count++;

                });

                if (!notSet)
                {
                    break;
                }
            }


            CHECK_EQUAL(solver->allPossibleReactions().at(choice), reaction);
            CHECK_EQUAL(choice, count);
            CHECK_EQUAL(choice, count2);

            count2++;

            r_pre = r_i;
        }

        n++;

    }
}



void testBed::testRateCalculation()
{

    double E, Esp;

    solver->setTargetConcentration(0.01);
    solver->initializeSolutionBath();

    solver->getRateVariables();

    for (SoluteParticle *particle : solver->particles())
    {
        particle->forEachActiveReactionDo([&] (Reaction * r)
        {

            if (!r->isType("DiffusionReaction"))
            {
                return;
            }

            E = ((DiffusionReaction*)r)->lastUsedEnergy();

            Esp = ((DiffusionReaction*)r)->lastUsedEsp();

            r->forceUpdateFlag(Reaction::defaultUpdateFlag);

            r->calcRate();

            CHECK_EQUAL(E, ((DiffusionReaction*)r)->lastUsedEnergy());

            CHECK_EQUAL(Esp, ((DiffusionReaction*)r)->lastUsedEsp());

        });
    }
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

        if (!currentSite->isActive())
        {
            return;
        }

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
            CHECK_EQUAL(nn(K), currentSite->associatedParticle()->nNeighbors(K));
        }

        CHECK_CLOSE(E, currentSite->associatedParticle()->energy(), 0.00001);

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
        CHECK_EQUAL(true, currentSite->isActive());

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
            CHECK_EQUAL(2*(12*(K+1)*(K+1) + 1), currentSite->associatedParticle()->nNeighbors(K) + nBlocked(K));
        }

        CHECK_CLOSE(eMax, currentSite->associatedParticle()->energy() + blockedE, 0.001);

        accuBlockedE += blockedE;

    });

    CHECK_EQUAL(NX()*NY()*NZ(), SoluteParticle::nParticles());
    CHECK_CLOSE(NX()*NY()*NZ()*eMax, SoluteParticle::totalEnergy() + accuBlockedE, 0.001);

    deactivateAllSites();

    CHECK_EQUAL(0, SoluteParticle::nParticles());
    CHECK_CLOSE(0, SoluteParticle::totalEnergy(), 0.001);

}


void testBed::testInitialReactionSetup()
{

    CHECK_EQUAL(0, SoluteParticle::nParticles());
    CHECK_EQUAL(0, Reaction::nReactions());

    solver->initializeCrystal(0.3);
    SoluteParticle::updateAffectedParticles();

    std::vector<Reaction*> oldReactions;
    double totRate1 = 0;

    for (SoluteParticle *particle : solver->particles())
    {
        particle->forEachActiveReactionDo([&] (Reaction * r)
        {
            oldReactions.push_back(r);
            totRate1 += r->rate();

        });
    }

    for (SoluteParticle *particle : solver->particles())
    {
        particle->forEachActiveReactionDo([] (Reaction * r)
        {
            r->forceUpdateFlag(Reaction::defaultUpdateFlag);
        });

        particle->updateReactions();

    }


    std::vector<Reaction*> reactions;
    double totRate2 = 0;

    for (SoluteParticle *particle : solver->particles())
    {
        particle->forEachActiveReactionDo([&] (Reaction * r)
        {
            reactions.push_back(r);
            totRate2 += r->rate();
        });
    }

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


    const SnapShot s00(solver);

    const SnapShot & s1 = *testSequentialCore();


    solver->reset();

    const SnapShot s01(solver);

    const SnapShot & s2 = *testSequentialCore();



    delete solver;

    makeSolver();

    initBoundaryTestParameters();

    const SnapShot s02(solver);

    const SnapShot & s3 = *testSequentialCore();


    delete solver;

    makeSolver();

    initBoundaryTestParameters();

    const SnapShot s03(solver);

    const SnapShot & s4 = *testSequentialCore();



    CHECK_EQUAL(s00, s01);
    CHECK_EQUAL(s01, s02);
    CHECK_EQUAL(s02, s03);
    CHECK_EQUAL(s01, s02);
    CHECK_EQUAL(s01, s03);
    CHECK_EQUAL(s02, s03);

    CHECK_EQUAL(s1, s2);
    CHECK_EQUAL(s3, s4); //!

    CHECK_EQUAL(s1, s3); //!
    CHECK_EQUAL(s2, s3); //!
    CHECK_EQUAL(s1, s4); //!
    CHECK_EQUAL(s2, s4); //!


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

    solver->setRNGSeed(Seed::specific, baseSeed);

    solver->setBoxSize({10, 10, 10}, false);

    Site::resetBoundariesTo(lastBoundaries);

    Site::resetNNeighborsLimitTo(3);

}

void testBed::initSimpleSystemParameters()
{

    solver->setBoxSize({15, 15, 15}, false);

    Site::resetNNeighborsLimitTo(2);

    Site::resetBoundariesTo(Boundary::Periodic);

}

void testBed::activateAllSites()
{

    solver->forEachSiteDo([] (Site * currentSite)
    {
        if (!currentSite->isActive())
        {
            SoluteParticle *particle = new SoluteParticle();
            solver->spawnParticle(particle, currentSite, false);
        }
    });

}

void testBed::deactivateAllSites()
{

    solver->forEachSiteDo([] (Site * currentSite)
    {
        if (currentSite->isActive())
        {
            solver->despawnParticle(currentSite);
        }
    });
}

Site *testBed::getBoxCenter(const int dx, const int dy, const int dz)
{
    return solver->getSite(NX()/2 + dx, NY()/2 + dy, NZ()/2 + dz);
}

void testBed::_reactionShufflerCheck(uint nReacs)
{

    CHECK_EQUAL(nReacs, solver->allPossibleReactions().size());
    CHECK_EQUAL(nReacs, solver->accuAllRates().size());
    CHECK_EQUAL(0,      solver->availableReactionSlots().size());

    for (uint i = 0; i < nReacs; ++i)
    {
        CHECK_EQUAL(i + 1, solver->accuAllRates().at(i));
        CHECK_EQUAL(i,     solver->allPossibleReactions().at(i)->address());
    }

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


void testBed::testOptimizedRateVectors()
{

    solver->setRNGSeed(Seed::specific, 1000000);

    Reaction::setLinearRateScale(1000);

    solver->initializeCrystal(0.2);

    double kTotBF;
    vector<double> accuAllRatesBF;
    vector<Reaction*> allPossibleReactionsBF;

    uint cycle = 1;
    uint nCycles = 500;

    bool containsReaction;

    while (cycle <= nCycles)
    {

        solver->getRateVariables();

        fill_rate_stuff(accuAllRatesBF, allPossibleReactionsBF, kTotBF);

        CHECK_CLOSE(kTotBF, solver->kTot(), 1E-5);
        CHECK_CLOSE(kTotBF, *(solver->accuAllRates().end() - 1), 1E-5);


        for (Reaction * r : solver->allPossibleReactions())
        {
            containsReaction = std::find(allPossibleReactionsBF.begin(), allPossibleReactionsBF.end(), r) != allPossibleReactionsBF.end();
            CHECK_EQUAL(true, containsReaction);
        }

        double prevAR = 0;
        for (double AR : solver->accuAllRates())
        {
            CHECK_EQUAL(true, prevAR < AR);
        }

        double kTotSBF = 0;
        for (uint i = 0; i < solver->allPossibleReactions().size(); ++i)
        {

            kTotSBF += solver->allPossibleReactions().at(i)->rate();

            CHECK_CLOSE(kTotSBF, solver->accuAllRates().at(i), 1E-5);

        }


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
    return;
    double c;

    Site * center = getBoxCenter();

    Reaction::setLinearRateScale(1000);


    //Activating the center, this should add 26 reactions with unit rate.
    solver->forceSpawnParticle(center);

    SoluteParticle::updateAffectedParticles();

    CHECK_EQUAL(0,  solver->availableReactionSlots().size());
    CHECK_EQUAL(26, solver->allPossibleReactions().size());
    CHECK_EQUAL(26, solver->accuAllRates().size());

    c = 0;
    for (uint i = 0; i < solver->allPossibleReactions().size(); ++i)
    {
        c += solver->allPossibleReactions().at(i)->rate();
        CHECK_EQUAL(c, solver->accuAllRates().at(i));
    }

    CHECK_EQUAL(solver->kTot(), *(solver->accuAllRates().end()-1));

    //Deactivating it should set all 26 reactions as available to overwrite.
    solver->despawnParticle(center);

    SoluteParticle::updateAffectedParticles();

    CHECK_EQUAL(26,  solver->availableReactionSlots().size());
    CHECK_EQUAL(26,  solver->allPossibleReactions().size());
    CHECK_EQUAL(26,  solver->accuAllRates().size());

    CHECK_EQUAL(0, solver->kTot());

    for (uint i = 0; i < solver->allPossibleReactions().size(); ++i)
    {
        CHECK_EQUAL(0, solver->accuAllRates().at(i));
    }

    //Activating again should reset back to original case. Not trivial because this
    //time the vacancy list is not empty.
    solver->forceSpawnParticle(center);

    SoluteParticle::updateAffectedParticles();

    CHECK_EQUAL(0,  solver->availableReactionSlots().size());
    CHECK_EQUAL(26, solver->allPossibleReactions().size());
    CHECK_EQUAL(26, solver->accuAllRates().size());

    CHECK_EQUAL(26*Reaction::linearRateScale(), solver->kTot());

    c = 0;
    for (uint i = 0; i < solver->allPossibleReactions().size(); ++i)
    {
        c += solver->allPossibleReactions().at(i)->rate();
        CHECK_EQUAL(c, solver->accuAllRates().at(i));
    }


    //activating a fixed crystal next to center particle. Fixed crystal has no
    //reactions so it should only induce 1 blocked reaction.
    Site * neighbor = getBoxCenter(1);
    solver->forceSpawnParticle(neighbor);

    SoluteParticle::updateAffectedParticles();

    CHECK_EQUAL(1,  solver->availableReactionSlots().size());

    c = 0;
    for (uint i = 0; i < solver->allPossibleReactions().size(); ++i)
    {
        if (i == solver->availableReactionSlots().at(0))
        {
            continue;
        }

        c += solver->allPossibleReactions().at(i)->rate();
        CHECK_CLOSE(c, solver->accuAllRates().at(i), 0.00001);
    }

    //activating a new particle independent of the others should now only induce 25 more spots,
    //since one is already vacant.

    Site * distantCousin = getBoxCenter(-2);

    SoluteParticle *p = new SoluteParticle();
    CHECK_EQUAL(true, solver->spawnParticle(p, distantCousin, true));

    SoluteParticle::updateAffectedParticles();

    CHECK_EQUAL(0,  solver->availableReactionSlots().size());
    CHECK_EQUAL(26+25, solver->allPossibleReactions().size());
    CHECK_EQUAL(26+25, solver->accuAllRates().size());

    c = 0;
    for (uint i = 0; i < solver->allPossibleReactions().size(); ++i)
    {
        c += solver->allPossibleReactions().at(i)->rate();
        CHECK_CLOSE(c, solver->accuAllRates().at(i), 0.00001);
    }


    //deactivating the neighbor should create two decoupled systems
    solver->despawnParticle(neighbor);

    SoluteParticle::updateAffectedParticles();

    CHECK_EQUAL(0,    solver->availableReactionSlots().size());
    CHECK_EQUAL(2*26, solver->allPossibleReactions().size());
    CHECK_EQUAL(2*26, solver->accuAllRates().size());

    CHECK_CLOSE(2*26*Reaction::linearRateScale(), solver->kTot(), 0.00001);

    c = 0;
    for (uint i = 0; i < solver->allPossibleReactions().size(); ++i)
    {
        c += solver->allPossibleReactions().at(i)->rate();
        CHECK_CLOSE(c, solver->accuAllRates().at(i), 0.00001);
    }

    //removing both the particles should bring us back to initial state

    solver->despawnParticle(center);
    solver->despawnParticle(distantCousin);

    SoluteParticle::updateAffectedParticles();

    CHECK_EQUAL(2*26, solver->availableReactionSlots().size());
    CHECK_EQUAL(2*26, solver->allPossibleReactions().size());
    CHECK_EQUAL(2*26, solver->accuAllRates().size());

    CHECK_CLOSE(0, solver->kTot(), 0.00001);

    for (uint i = 0; i < solver->allPossibleReactions().size(); ++i)
    {
        CHECK_CLOSE(0, solver->accuAllRates().at(i), 0.00001);
    }


    Reaction::setLinearRateScale(1);

}

void testBed::testReactionShuffler()
{
    return;
    uint nReacs = 100;

    vector<DummyReaction*> allReacs;

    //Initialize nReacs reactions and check that accuallrates is equal to 1, 2, 3 ...
    //and that the addresses are correct.
    for (uint i = 0; i < nReacs; ++i)
    {
        allReacs.push_back(new DummyReaction(i));
    }


    for (Reaction * r : allReacs)
    {
        r->calcRate();
    }

    SoluteParticle::clearAffectedParticles();

    _reactionShufflerCheck(nReacs);

    //this should do nothing: Everything is ordered.
    solver->reshuffleReactions();

    _reactionShufflerCheck(nReacs);

    allReacs.at(nReacs/2)->allowed = false;
    allReacs.at(nReacs/2)->disable();

    //this should replace the disabled reaction with the end reaction: Everything is ordered.
    solver->reshuffleReactions();

    CHECK_EQUAL(nReacs - 1, static_cast<DummyReaction*>(solver->allPossibleReactions().at(nReacs/2))->initAddress);

    _reactionShufflerCheck(nReacs - 1);

    //resetting, so far we have swapped the last and the middle element.
    SoluteParticle::clearAffectedParticles();

    allReacs.at(nReacs - 1)->allowed = false;
    allReacs.at(nReacs - 1)->disable();

    allReacs.at(nReacs/2)->allowed = true;
    allReacs.at(nReacs/2)->calcRate();

    SoluteParticle::clearAffectedParticles();

    allReacs.at(nReacs - 1)->allowed = true;
    allReacs.at(nReacs - 1)->calcRate();

    SoluteParticle::clearAffectedParticles();

    _reactionShufflerCheck(nReacs);

    for (DummyReaction * r : allReacs)
    {
        CHECK_EQUAL(r->address(), r->initAddress);
    }

    //remove every other element

    for (uint i = 0; i < nReacs; i += 2)
    {
        allReacs.at(i)->allowed = false;
        allReacs.at(i)->disable();
    }

    CHECK_EQUAL(nReacs/2, solver->availableReactionSlots().size());

    solver->reshuffleReactions();

    _reactionShufflerCheck(nReacs/2);


    //remove the rest of the elements
    for (uint i = 1; i < nReacs; i += 2)
    {
        allReacs.at(i)->allowed = false;
        allReacs.at(i)->disable();
    }

    CHECK_EQUAL(nReacs/2, solver->availableReactionSlots().size());

    for (uint addr : solver->availableReactionSlots())
    {
        CHECK_EQUAL(Reaction::UNSET_ADDRESS, solver->allPossibleReactions().at(addr)->address());
    }

    solver->reshuffleReactions();

    _reactionShufflerCheck(0);

    SoluteParticle::clearAffectedParticles();

    for (uint i = 0; i < nReacs; ++i)
    {
        allReacs.at(i)->allowed = true;
        allReacs.at(i)->calcRate();
    }

    CHECK_EQUAL(0, solver->availableReactionSlots().size());

    solver->reshuffleReactions();

    _reactionShufflerCheck(nReacs);

    for (DummyReaction * r : allReacs)
    {
        CHECK_EQUAL(r->address(), r->initAddress);
    }

    SoluteParticle::clearAffectedParticles();


    //Something is active, becomes deactivated, but is initiated to take
    //up a vacant spot before it is itself vacated?
    //YES SHUFFLE HAPPENS BEFORE VACATING NEW ONES? NO. Why is the address then not in vacant? ADDRESS IS SET WRONGLY.. what induses this?

    bool remove = false;
    uint c = 0;
    for (uint i = 0; i < nReacs; ++i)
    {

        if (remove)
        {
            allReacs.at(i)->allowed = false;
            allReacs.at(i)->disable();
            c++;
        }

        if (i % 5 == 0)
        {
            remove = !remove;
        }

    }

    CHECK_EQUAL(c, solver->availableReactionSlots().size());

    solver->reshuffleReactions();

    SoluteParticle::clearAffectedParticles();

    _reactionShufflerCheck(nReacs - c);


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

    for (SoluteParticle *particle : solver->particles())
    {
        particle->forEachActiveReactionDo([&] (Reaction * reaction)
        {

            KMCDebugger_Assert(reaction->rate(), !=, Reaction::UNSET_RATE, "Reaction rate should not be unset at this point.", reaction->getFinalizingDebugMessage());

            kTot += reaction->rate();

            accuAllRates.push_back(kTot);

            allPossibleReactions.push_back(reaction);

        });
    }
}


void testBed::testStateChanges()
{
    Site * center = getBoxCenter();

    solver->forceSpawnParticle(center);

    CHECK_EQUAL(ParticleStates::solvant, center->associatedParticle()->particleState());

    solver->despawnParticle(center);

    for (int i = -1; i <= 1; ++i)
    {
        for (int j = -1; j <= 1; ++j)
        {
            for (int k = -1; k <= 1; ++k)
            {
                solver->forceSpawnParticle(getBoxCenter(i, j, k));
            }
        }
    }

    solver->forceSpawnParticle(getBoxCenter(3, 0, 0));

    for (int i = -1; i <= 1; ++i)
    {
        for (int j = -1; j <= 1; ++j)
        {
            for (int k = -1; k <= 1; ++k)
            {
                if (i == j && j == k && k == 0)
                {
                    CHECK_EQUAL(ParticleStates::crystal, getBoxCenter(i, j, k)->associatedParticle()->particleState());
                }

                else
                {
                    CHECK_EQUAL(ParticleStates::surface, getBoxCenter(i, j, k)->associatedParticle()->particleState());
                }
            }
        }
    }

    CHECK_EQUAL(ParticleStates::solvant, getBoxCenter(3, 0, 0)->associatedParticle()->particleState());

}

void testBed::testNeighborlist()
{
    solver->setBoxSize({15, 15, 15});

    Site::resetNNeighborsLimitTo(5);

    solver->forceSpawnParticle(getBoxCenter());

    auto NNSUMBF = [&] (Site * site) -> uint
    {
        uint nn = 0;
        for (uint i = 0; i < Site::nNeighborsLimit(); ++i)
        {
            nn += site->associatedParticle()->neighbouringParticles().at(i).size();
        }

        return nn;
    };

    CHECK_EQUAL(0, NNSUMBF(getBoxCenter()));

    uint c = 0;
    uvec nn(Site::nNeighborsLimit(), fill::zeros);

    getBoxCenter()->forEachNeighborDo_sendIndices([&] (Site *site, uint i, uint j, uint k)
    {
        uint level = Site::levelMatrix(i, j, k);

        solver->forceSpawnParticle(site);

        c++;
        nn(level)++;

        CHECK_EQUAL(getBoxCenter()->associatedParticle()->isNeighbor(site->associatedParticle(), level), true);
        CHECK_EQUAL(site->associatedParticle()->isNeighbor(getBoxCenter()->associatedParticle(), level), true);

        CHECK_EQUAL(c, NNSUMBF(getBoxCenter()));
        CHECK_EQUAL(nn(level), getBoxCenter()->associatedParticle()->neighbouringParticles(level).size());

    });

    getBoxCenter()->forEachNeighborDo_sendIndices([&] (Site *site, uint i, uint j, uint k)
    {
        uint level = Site::levelMatrix(i, j, k);

        solver->despawnParticle(site);

        c--;
        nn(level)--;

        CHECK_EQUAL(getBoxCenter()->associatedParticle()->isNeighbor(site->associatedParticle(), level), false);


        CHECK_EQUAL(c, NNSUMBF(getBoxCenter()));
        CHECK_EQUAL(nn(level), getBoxCenter()->associatedParticle()->neighbouringParticles(level).size());

    });

}

void testBed::testAccuAllRates()
{

    solver->initializeCrystal(0.3);

    solver->getRateVariables();

    uint i = 0;
    double kTot = 0;

    for (Reaction *r : solver->allPossibleReactions())
    {
        kTot += r->rate();

        CHECK_CLOSE(kTot, solver->accuAllRates().at(i), solver->minRateThreshold());

        i++;

    }


}

void testBed::testInitialSiteSetup()
{

    int nx, ny, nz;

    uint xt, yt, zt;

    Site *site;



    CHECK_EQUAL(NX()*NY()*NZ(), Site::_refCount());

    solver->clearSites();

    CHECK_EQUAL(0, Site::_refCount());

    solver->initializeSites();

    CHECK_EQUAL(NX()*NY()*NZ(), Site::_refCount());


    solver->forEachSiteDo([&] (Site *_site)
    {
       site = _site;

       CHECK_EQUAL(site, solver->getSite(site->x(), site->y(), site->z()));

       CHECK_EQUAL(site->x(), solver->getSite(site->x(), site->y(), site->z())->x());
       CHECK_EQUAL(site->y(), solver->getSite(site->x(), site->y(), site->z())->y());
       CHECK_EQUAL(site->z(), solver->getSite(site->x(), site->y(), site->z())->z());
    });

    Site::resetBoundariesTo(Boundary::Periodic);


    solver->forEachSiteDo([&] (Site *_site)
    {
        site = _site; //hack to fix autocompletion.

        Boundary::setupCurrentBoundaries(site->x(), site->y(), site->z());

        for (uint i = 0; i < Site::neighborhoodLength(); ++i)
        {
            nx = (int)site->x() + Site::originTransformVector(i);
            xt = Boundary::currentBoundaries(0)->transformCoordinate(nx);

            for (uint j = 0; j < Site::neighborhoodLength(); ++j)
            {
                ny = (int)site->y() + Site::originTransformVector(j);
                yt = Boundary::currentBoundaries(1)->transformCoordinate(ny);

                for (uint k = 0; k < Site::neighborhoodLength(); ++k)
                {
                    nz = (int)site->z() + Site::originTransformVector(k);
                    zt = Boundary::currentBoundaries(2)->transformCoordinate(nz);

                    CHECK_EQUAL(false, site->neighborhood(i, j, k) == NULL);
                    CHECK_EQUAL(site->neighborhood(i, j, k), solver->getSite(nx, ny, nz));

                    if ((nx < 0 || nx >= (int)NX()) ||
                        (ny < 0 || ny >= (int)NY()) ||
                        (nz < 0 || nz >= (int)NZ()))
                    {
                        CHECK_EQUAL(*solver->getSite(nx, ny, nz), *solver->getSite(xt, yt, zt));
                        cout.flush();
                    }

                }
            }
        }

    });

    Site::resetBoundariesTo(Boundary::Edge);

    solver->forEachSiteDo([&] (Site *_site)
    {
        site = _site; //hack to fix autocompletion.

        Boundary::setupCurrentBoundaries(site->x(), site->y(), site->z());

        for (uint i = 0; i < Site::neighborhoodLength(); ++i)
        {
            nx = (int)site->x() + Site::originTransformVector(i);
            xt = Boundary::currentBoundaries(0)->transformCoordinate(nx);

            for (uint j = 0; j < Site::neighborhoodLength(); ++j)
            {
                ny = (int)site->y() + Site::originTransformVector(j);
                yt = Boundary::currentBoundaries(1)->transformCoordinate(ny);

                for (uint k = 0; k < Site::neighborhoodLength(); ++k)
                {
                    nz = (int)site->z() + Site::originTransformVector(k);
                    zt = Boundary::currentBoundaries(2)->transformCoordinate(nz);

                    CHECK_EQUAL(site->neighborhood(i, j, k), solver->getSite(nx, ny, nz));

                    if (solver->getSite(nx, ny, nz) == NULL)
                    {
                        CHECK_EQUAL(true, Boundary::isBlocked(xt, yt, zt));
                    }

                }
            }
        }

    });


}


KMCSolver * testBed::solver;

wall_clock testBed::timer;

seed_type testBed::baseSeed;


string testBed::lastBoundariesName;

umat testBed::lastBoundaries;
