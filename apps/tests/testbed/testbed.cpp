#include "testbed.h"

#include "../snapshot/snapshot.h"

#include "../dummyreaction.h"

#include <unittest++/UnitTest++.h>

#include <iostream>

void testBed::makeSolver()
{

    solver = new KMCSolver();

    initSimpleSystemParameters(false);

}

void testBed::testTotalParticleStateCounters()
{

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

    solver->reset();

    CHECK_EQUAL(0, accu(SoluteParticle::totalParticlesVector()));

    uint C1 = 0;

    solver->forEachSiteDo([&C1] (uint x, uint y, uint z, Site * currentSite)
    {
        (void) currentSite;

        CHECK_EQUAL(C1, accu(SoluteParticle::totalParticlesVector()));

        solver->forceSpawnParticle(x, y, z);
        C1++;
    });

    solver->forEachSiteDo([&C1] (uint x, uint y, uint z, Site * currentSite)
    {
        (void) x;
        (void) y;
        (void) z;

        C1--;

        solver->despawnParticle(currentSite->associatedParticle());

        CHECK_EQUAL(C1, accu(SoluteParticle::totalParticlesVector()));

    });

    CHECK_EQUAL(0, accu(SoluteParticle::totalParticlesVector()));

    forceNewBoxSize({10, 10, 10});

    solver->forEachSiteDo([] (uint x, uint y, uint z, Site * site)
    {
        CHECK_EQUAL(false, site->isActive());

        if ((x%2 == 0) && (y%2 == 0) && (z%2 == 0))
        {
            solver->forceSpawnParticle(x, y, z);
        }

    });

    CHECK_EQUAL(NX()*NY()*NZ()/8, SoluteParticle::totalParticles(ParticleStates::solvant));



}

void testBed::testDistanceTo()
{

    int dx, dy, dz, dx2, dy2, dz2;
    uint adx, ady, adz;

    forceNewBoxSize({8, 8, 8});

    solver->forEachSiteDo([&] (uint startx, uint starty, uint startz, Site * startSite)
    {
        (void) startSite;

        solver->forEachSiteDo([&] (uint endx, uint endy, uint endz, Site * endSite)
        {

            (void) endSite;

            Boundary::setupCurrentBoundaries(endx, endy, endz);

            Site::distanceBetween(endx, endy, endz, startx, starty, startz, dx, dy, dz, true);

            adx = dx;
            ady = dy;
            adz = dz;

            Site::distanceBetween(endx, endy, endz, startx, starty, startz, dx, dy, dz);

            CHECK_EQUAL(adx, abs(dx));
            CHECK_EQUAL(ady, abs(dy));
            CHECK_EQUAL(adz, abs(dz));

            Site::distanceBetween(startx, starty, startz, endx, endy, endz, dx2, dy2, dz2);

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

            CHECK_EQUAL(startx, Boundary::currentBoundaries(0)->transformCoordinate(endx + dx));
            CHECK_EQUAL(starty, Boundary::currentBoundaries(1)->transformCoordinate(endy + dy));
            CHECK_EQUAL(startz, Boundary::currentBoundaries(2)->transformCoordinate(endz + dz));

            if (Site::boundaryTypes(0) == Boundary::Periodic)
            {
                int X = (int)startx - (int)endx;

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
                int Y = (int)starty - (int)endy;

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
                int Z = (int)startz - (int)endz;

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

                    const Site & site = *(currentDiffReaction->reactant()->site());
                    const Site & dest  = *(currentDiffReaction->destinationSite());

                    CHECK_EQUAL(particle->site(), &site);

                    Site::distanceBetween(site.associatedParticle()->x(),
                                          site.associatedParticle()->y(),
                                          site.associatedParticle()->z(),
                                          site.associatedParticle()->x() + currentDiffReaction->path(0),
                                          site.associatedParticle()->y() + currentDiffReaction->path(1),
                                          site.associatedParticle()->z() + currentDiffReaction->path(2),
                                          _i, _j, _k);

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

    solver->forEachSiteDo([&] (uint x, uint y, uint z, Site * currentSite)
    {
        (void) currentSite;

        Site::forEachNeighborDo_sendPath(x, y, z, [&] (Site * neighbor, int i, int j, int k)
        {
            (void) neighbor;

            Site::distanceBetween(x, y, z, x + i, y + j, z + k, dx, dy, dz);

            CHECK_EQUAL(i, dx);
            CHECK_EQUAL(j, dy);
            CHECK_EQUAL(k, dz);

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

    solver->forceSpawnParticle(0, 0, 0);
    solver->forceSpawnParticle(2, 0, 0);
    solver->forceSpawnParticle(0, 2, 0);


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

    vector<vector<uint> >crystalSites = {{sx, sy    , sz    },

                                         {sx, sy    , sz + 1},
                                         {sx, sy    , sz + 2},

                                         {sx, sy - 1, sz - 1},
                                         {sx, sy - 1, sz    },
                                         {sx, sy - 1, sz + 1},

                                         {sx, sy + 1, sz + 1}
                                        };

    Site *currentSite;


    for (auto xyz: crystalSites)
    {
        solver->forceSpawnParticle(xyz[0], xyz[1], xyz[2]);
    }

    for (auto xyz: crystalSites)
    {

        for (int i = -1; i <= 1; ++i)
        {
            for (int j = -1; j <= 1; ++j)
            {
                for (int k = -1; k <= 1; ++k)
                {

                    currentSite = solver->getSite(i + xyz[0], j + xyz[1], k + xyz[2]);

                    if (!currentSite->isActive())
                    {
                        solver->forceSpawnParticle(i + xyz[0], j + xyz[1], k + xyz[2]);
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

void testBed::testParticleMixing()
{

    SoluteParticle *A, *B;
    DiffusionReaction *rA, *rB;

    forceNewBoundaries(Boundary::Edge); //Not periodic
    forceNewBoxSize({10, 1, 1});        //make a 1D system

    vector<double> strengths = {1.0, 10.0, 100.0};
    vector<double> powers = {1.0, 2.0, 3.0};

    DiffusionReaction::setPotentialParameters(powers, strengths);

    CHECK_EQUAL(SoluteParticle::nSpecies(), strengths.size());

    for (uint typeA = 0; typeA < SoluteParticle::nSpecies(); ++typeA)
    {
        //Spawning a particle of type A. Since it's a 1D system, it should have 2 reactions.
        A = forceSpawnCenter(-1, 0, 0, typeA);

        CHECK_EQUAL(2, A->reactions().size());
        rA = A->diffusionReactions(2, 0, 0); //left-going reaction

        for (uint typeB = 0; typeB < SoluteParticle::nSpecies(); ++typeB)
        {
            B = forceSpawnCenter(1, 0, 0, typeB); //A and B are at a distance 2 from eachother.

            solver->getRateVariables();

            CHECK_EQUAL(2, B->reactions().size());
            rB = B->diffusionReactions(0, 0, 0); //right->going reaction

            CHECK_EQUAL(A->energy(), B->energy());

            //checking value of energy
            double rPowerCombo = sqrt(powers.at(typeA)*powers.at(typeB));
            double strengthCombo = 0.5*(strengths.at(typeA) + strengths.at(typeB));
            double eCombo = strengthCombo/pow(2.0, rPowerCombo);

            CHECK_CLOSE(eCombo, A->energy(), 1E-10);

            //and the saddle energies of the reactions of moving towards the center (site between A and B)
            CHECK_EQUAL(rA->lastUsedEsp(), rB->lastUsedEsp());


            solver->despawnParticle(B);

        }

        solver->despawnParticle(A);
    }



}

void testBed::testBoundarySites()
{
    cout << "test outdated." << endl;
    return;

    forceNewBoxSize({10, 15, 20});

    uvec sizes = {NY()*NZ(), NX()*NZ(), NX()*NY()};

    vector<vector<uint> > shapes = {{NY(), NZ()}, {NX(), NZ()}, {NX(), NY()}};

    vector<vector<uint> > indices = {{1, 2}, {0, 2}, {0, 1}};

    vector<uint> shape;

    uvec xyz;

    uint size;

    uint x, y, z;

    uint c1, c2;

    for (uint dim = 0; dim < 3; ++dim)
    {
        size = sizes(dim);

        shape = shapes.at(dim);


        for (uint orientation = 0; orientation < 2; ++orientation)
        {
            CHECK_EQUAL(size, Site::boundaries(dim, orientation)->boundarySize());

            c1 = 0;
            c2 = 0;

            for (uint n = 0; n < size; ++n)
            {
                Site::boundaries(dim, orientation)->getBoundarySite(n, x, y, z);

                xyz = {x, y, z};

                if (c1 == shape.at(1))
                {
                    c1 = 0;
                    c2++;
                }

                CHECK_EQUAL(Site::boundaries(dim, orientation)->bound(), xyz(dim));

                CHECK_EQUAL(c1, xyz(indices.at(dim).at(1)));
                CHECK_EQUAL(c2, xyz(indices.at(dim).at(0)));

                c1++;
            }

        }
    }
}

void testBed::testConcentrationWall()
{
    cout << "test outdated." << endl;
    return;

    forceNewBoxSize({10, 10, 10});

    uint outerShellSize = 5*5*5 - 3*3*3;
    uint size = 10*10/4;

    forceNewBoundaries(Boundary::ConcentrationWall);

    forceNewNNeighborLimit(1);

    solver->setTargetConcentration(outerShellSize/((double)NX()*NY()*NZ()));

    ConcentrationWall::setMaxEventsPrCycle(100);

    CHECK_EQUAL(0, SoluteParticle::nParticles());


    Site::updateBoundaries();

    CHECK_EQUAL(outerShellSize, SoluteParticle::nParticles());

    ConcentrationWall::setMaxEventsPrCycle(size/2);

    solver->setTargetConcentration(solver->targetConcentration()/2);

    Site::updateBoundaries();

    CHECK_EQUAL(outerShellSize/2, SoluteParticle::nParticles());



    ConcentrationWall::setMaxEventsPrCycle(1);

    solver->setTargetConcentration(solver->targetConcentration()*2);

    Site::updateBoundaries();

    CHECK_EQUAL(outerShellSize/2 + 6, SoluteParticle::nParticles());


    ConcentrationWall::setMaxEventsPrCycle(100);

    solver->setTargetConcentration(0);

    Site::updateBoundaries();

    CHECK_EQUAL(0, SoluteParticle::nParticles());


    ConcentrationWall::setMaxEventsPrCycle(3);
    solver->setTargetConcentration(0.01);

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

#ifdef NDEBUG
    cout << "test requires debug mode. ";
    return;
#endif


    uint choice, count, count2;

    double kTot, r_pre;

    Reaction* reaction;


    uint N = 2;
    uint n = 0;

    KMCDebugger_Assert(0, ==, SoluteParticle::affectedParticles().size());
    KMCDebugger_Assert(0, ==, solver->allPossibleReactions().size());
    KMCDebugger_Assert(0, ==, solver->accuAllRates().size());

    solver->setTargetConcentration(10./(NX()*NY()*NZ()));
    solver->initializeCrystal(0.2);

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

                    CHECK_EQUAL(true, solver->isPossibleReaction(r));
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

void testBed::testAffectedParticles()
{
    solver->initializeSolutionBath();

    for (SoluteParticle *particle : solver->particles())
    {
        CHECK_EQUAL(true, particle->isAffected());
    }

    solver->getRateVariables();

    for (SoluteParticle *particle : solver->particles())
    {
        CHECK_EQUAL(false, particle->isAffected());
    }


    SoluteParticle *particle = solver->particle(0);

    for (Reaction *r : particle->reactions())
    {
        if (r->isAllowed())
        {
            r->execute();
            break;
        }
    }

    CHECK_EQUAL(true, particle->isAffected());

    solver->despawnParticle(particle);

    CHECK_EQUAL(false, SoluteParticle::isAffected(particle));

    deactivateAllSites();

    Site *center = getBoxCenter();
    forceSpawnCenter();

    CHECK_EQUAL(true, center->associatedParticle()->isAffected());
    CHECK_EQUAL(1, SoluteParticle::affectedParticles().size());

    solver->getRateVariables();

    CHECK_EQUAL(false, center->associatedParticle()->isAffected());
    CHECK_EQUAL(0, SoluteParticle::affectedParticles().size());

    Site *neighbor = getBoxCenter(1);
    forceSpawnCenter(1);

    CHECK_EQUAL(true, center->associatedParticle()->isAffected());
    CHECK_EQUAL(true, neighbor->associatedParticle()->isAffected());
    CHECK_EQUAL(2, SoluteParticle::affectedParticles().size());

    solver->getRateVariables();

    CHECK_EQUAL(false, center->associatedParticle()->isAffected());
    CHECK_EQUAL(false, neighbor->associatedParticle()->isAffected());
    CHECK_EQUAL(0, SoluteParticle::affectedParticles().size());

    center->associatedParticle()->markAsAffected();

    CHECK_EQUAL(true, center->associatedParticle()->isAffected());
    CHECK_EQUAL(1, SoluteParticle::affectedParticles().size());

    center->associatedParticle()->markAsAffected();
    center->associatedParticle()->markAsAffected();
    center->associatedParticle()->markAsAffected();
    center->associatedParticle()->markAsAffected();
    center->associatedParticle()->markAsAffected();

    CHECK_EQUAL(true, center->associatedParticle()->isAffected());
    CHECK_EQUAL(1, SoluteParticle::affectedParticles().size());

    for (Reaction *r : center->associatedParticle()->reactions())
    {
        r->forceUpdateFlag(Reaction::defaultUpdateFlag);
    }

    solver->getRateVariables();

    CHECK_EQUAL(false, center->associatedParticle()->isAffected());
    CHECK_EQUAL(0, SoluteParticle::affectedParticles().size());


}



void testBed::testRateCalculation()
{

    double E, Esp;

    forceNewNNeighborLimit(2);
    solver->setTargetConcentration(0.01);
    solver->initializeSolutionBath();

    solver->getRateVariables();

    DiffusionReaction *r;

    for (SoluteParticle *particle : solver->particles())
    {
        particle->forEachActiveReactionDo([&] (Reaction * _r)
        {

            if (!_r->isType("DiffusionReaction"))
            {
                return;
            }

            r = (DiffusionReaction*)_r;

            E = r->lastUsedEnergy();

            Esp = r->lastUsedEsp();

            r->forceUpdateFlag(Reaction::defaultUpdateFlag);

            r->calcRate();

            CHECK_EQUAL(E, r->lastUsedEnergy());

            CHECK_EQUAL(Esp, r->lastUsedEsp());

        });
    }
}

void testBed::testEnergyAndNeighborSetup()
{

    int dx, dy, dz;
    uint ldx, ldy, ldz;

    double E;
    uint C;

    initSimpleSystemParameters();

    solver->initializeCrystal(0.3);

    uvec nn(Site::nNeighborsLimit());


    solver->forEachSiteDo([&] (uint x0, uint y0, uint z0, Site * currentSite)
    {

        if (!currentSite->isActive())
        {
            return;
        }

        E = 0;
        C = 0;
        nn.zeros();

        solver->forEachSiteDo([&] (uint x1, uint y1, uint z1, Site * otherSite)
        {

            Site::distanceBetween(x0, y0, z0, x1, y1, z1, dx, dy, dz);

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

                        CHECK_EQUAL(otherSite, Site::neighborhood(x0, y0, z0, dx, dy, dz));

                        CHECK_EQUAL(currentSite, Site::neighborhood(x1, y1, z1, -dx, -dy, -dz));
                    }
                }
            }
        });

        uint nNeighbors = 0;
        Site::forEachNeighborDo(x0, y0, z0, [&nNeighbors] (Site * neighbor)
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

    double eMax = 0;
    Site::forEachNeighborDo_sendIndices(0, 0, 0, [&eMax] (Site * site, uint i, uint j, uint k)
    {
        (void) site;
        eMax += DiffusionReaction::potential(i, j, k);
    });

    double blockedE;
    uvec nBlocked(Site::nNeighborsLimit());

    double accuBlockedE = 0;

    solver->forEachSiteDo([&] (uint x0, uint y0, uint z0, Site * currentSite)
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
                    if (Site::neighborhood_fromIndex(x0, y0, z0, ni, nj, nk) == NULL)
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

#ifdef NDEBUG
    cout << "test requires debug mode. ";
    return;
#endif

    solver->reset();

    initBoundaryTestParameters();


    const SnapShot s00(solver);

    const SnapShot & s1 = *sequentialCore();


    solver->reset();

    const SnapShot s01(solver);

    const SnapShot & s2 = *sequentialCore();



    delete solver;

    makeSolver();

    initBoundaryTestParameters();

    const SnapShot s02(solver);

    const SnapShot & s3 = *sequentialCore();


    delete solver;

    makeSolver();

    initBoundaryTestParameters();

    const SnapShot s03(solver);

    const SnapShot & s4 = *sequentialCore();



    CHECK_EQUAL(s00, s01);
    CHECK_EQUAL(s01, s02);
    CHECK_EQUAL(s02, s03);
    CHECK_EQUAL(s01, s02);
    CHECK_EQUAL(s01, s03);
    CHECK_EQUAL(s02, s03);

    CHECK_EQUAL(s1, s2);
    CHECK_EQUAL(s3, s4);

    CHECK_EQUAL(s1, s3); //!
    CHECK_EQUAL(s2, s3); //!
    CHECK_EQUAL(s1, s4); //!
    CHECK_EQUAL(s2, s4); //!


}

const SnapShot * testBed::sequentialCore()
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


    solver->clearSites();


    solver->setBoxSize({10, 10, 10}, false);

    Site::resetBoundariesTo(lastBoundaries);

    Site::resetNNeighborsLimitTo(3);


    solver->initializeSites();

    Site::initializeBoundaries();

}

void testBed::initSimpleSystemParameters(bool clean)
{

    if (clean)
    {
        solver->clearSites();
    }

    solver->setRNGSeed();


    solver->setNumberOfCycles(1000);

    solver->setTargetConcentration(0.005);

    Reaction::setBeta(0.5);

    DiffusionReaction::setPotentialParameters({1.0}, {1.0}, false);


    solver->setBoxSize({15, 15, 15}, false);

    if (clean)
    {
        Site::resetNNeighborsLimitTo(2);

        Site::resetBoundariesTo(Boundary::Periodic);
    }

    else
    {
        Site::setInitialNNeighborsLimit(2);

        Site::setInitialBoundaries(Boundary::Periodic);
    }


    solver->initializeSites();

    Site::initializeBoundaries();

}

void testBed::activateAllSites()
{

    solver->forEachSiteDo([] (uint x, uint y, uint z, Site * currentSite)
    {

        if (!currentSite->isActive())
        {
            SoluteParticle *particle = new SoluteParticle();
            solver->spawnParticle(particle, x, y, z, false);
        }
    });

}

void testBed::deactivateAllSites()
{

    solver->forEachSiteDo([] (uint x, uint y, uint z, Site * currentSite)
    {
        (void) x;
        (void) y;
        (void) z;

        if (currentSite->isActive())
        {
            solver->despawnParticle(currentSite->associatedParticle());
        }
    });
}

Site *testBed::getBoxCenter(const int dx, const int dy, const int dz)
{
    return solver->getSite(NX()/2 + dx, NY()/2 + dy, NZ()/2 + dz);
}

SoluteParticle *testBed::forceSpawnCenter(const int dx, const int dy, const int dz, const uint particleType)
{
    return solver->forceSpawnParticle(NX()/2 + dx, NY()/2 + dy, NZ()/2 + dz, particleType);
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

void testBed::forceNewBoxSize(const uvec3 boxSize, bool check)
{
    Site::finalizeBoundaries();

    solver->clearSites();
    solver->setBoxSize(boxSize, check);
    solver->initializeSites();

    Site::initializeBoundaries();
}

void testBed::forceNewNNeighborLimit(const uint nNeighborlimit, bool check)
{
    solver->clearSites();
    Site::resetNNeighborsLimitTo(nNeighborlimit, check);
    solver->initializeSites();
}

void testBed::forceNewBoundaries(const umat &boundaryMatrix)
{
    solver->clearSites();
    Site::resetBoundariesTo(boundaryMatrix);
    solver->initializeSites();

    Site::initializeBoundaries();
}

void testBed::forceNewBoundaries(const int boundaryType)
{
    forceNewBoundaries(umat::fixed<3, 2>(fill::zeros) + boundaryType);
}

void testBed::testKnownCase()
{

#ifdef NDEBUG
    cout << "test requires debug mode. ";
    return;
#endif

    delete solver;

    Config cfg;

    cfg.readFile("infiles/knowncase.cfg");

    const Setting & root = cfg.getRoot();

    solver = new KMCSolver(root);

    forceNewBoundaries(lastBoundaries);

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

    solver->forEachSiteDo([&] (uint x, uint y, uint z, Site * site)
    {
        (void) x;
        (void) y;
        (void) z;

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

    forceNewNNeighborLimit(2, false);

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

                forceNewBoxSize(boxSize);

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

                solver->forEachSiteDo([&] (uint x0, uint y0, uint z0, Site * currentSite)
                {
                    (void) currentSite;

                    nBlocked = 0;

                    for (const int &dx : Site::originTransformVector())
                    {
                        for (const int &dy : Site::originTransformVector())
                        {
                            for (const int &dz : Site::originTransformVector())
                            {
                                if (Site::neighborhood(x0, y0, z0, dx, dy, dz) == NULL)
                                {
                                    nBlocked++;
                                }
                            }
                        }
                    }


                    uint nNeighbors = 0;
                    Site::forEachNeighborDo(x0, y0, z0, [&nNeighbors] (Site * neighbor)
                    {
                        (void) neighbor;
                        nNeighbors++;
                    });

                    CHECK_EQUAL(pow(Site::neighborhoodLength(), 3) - 1, nNeighbors + nBlocked);

                    Site::forEachNeighborDo(x0, y0, z0, [&allSites] (Site * neighbor)
                    {
                        allSites.insert(neighbor);
                    });

                });

                CHECK_EQUAL(nx*ny*nz, allSites.size());

            }
        }
    }

}

void testBed::testnNeighborsLimit()
{


    set<Site*> allSites;

    uvec nNlims = {1, 2, 3};

    uvec3 boxSize = {10, 10, 10};


    solver->clearSites();

    Site::resetBoundariesTo(Boundary::Periodic);

    solver->setBoxSize(boxSize, false);

    solver->initializeSites();

    Site::initializeBoundaries();


    for (uint nNlim : nNlims)
    {
        allSites.clear();
        forceNewNNeighborLimit(nNlim, false);

        solver->forEachSiteDo([&] (uint x0, uint y0, uint z0, Site * currentSite)
        {
            (void) currentSite;

            Site::forEachNeighborDo(x0, y0, z0, [&allSites] (Site * neighbor)
            {
                allSites.insert(neighbor);
            });

            uint nNeighbors = 0;
            Site::forEachNeighborDo(x0, y0, z0, [&nNeighbors] (Site * neighbor)
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

#ifdef NDEBUG
    cout << "test requires debug mode. ";
    return;
#endif

    double c;

    Site * center = getBoxCenter();

    Reaction::setLinearRateScale(1000);


    //Activating the center, this should add 26 reactions with unit rate.
    forceSpawnCenter();

    SoluteParticle::updateAffectedParticles();

    CHECK_EQUAL(0,  solver->availableReactionSlots().size());
    CHECK_EQUAL(26, solver->allPossibleReactions().size());
    CHECK_EQUAL(26, solver->accuAllRates().size());

    c = 0;
    for (uint i = 0; i < solver->allPossibleReactions().size(); ++i)
    {
        CHECK_EQUAL(Reaction::linearRateScale(), solver->allPossibleReactions().at(i)->rate());

        c += solver->allPossibleReactions().at(i)->rate();
        CHECK_EQUAL(c, solver->accuAllRates().at(i));
    }

    CHECK_EQUAL(solver->kTot(), *(solver->accuAllRates().end()-1));

    //Deactivating it should set all 26 reactions as available to overwrite. Removes it from affected particles.
    solver->despawnParticle(center->associatedParticle());

    CHECK_EQUAL(26,  solver->availableReactionSlots().size());
    CHECK_EQUAL(26,  solver->allPossibleReactions().size());
    CHECK_EQUAL(26,  solver->accuAllRates().size());

    CHECK_EQUAL(0, solver->kTot());

    for (uint i = 0; i < solver->allPossibleReactions().size(); ++i)
    {
        CHECK_EQUAL(true, solver->isEmptyAddress(i));
        CHECK_EQUAL(0, solver->accuAllRates().at(i));
    }

    //Activating again should reset back to original case. Not trivial because this
    //time the vacancy list is not empty.
    forceSpawnCenter();

    SoluteParticle::updateAffectedParticles();

    CHECK_EQUAL(0,  solver->availableReactionSlots().size());
    CHECK_EQUAL(26, solver->allPossibleReactions().size());
    CHECK_EQUAL(26, solver->accuAllRates().size());

    CHECK_EQUAL(26*Reaction::linearRateScale(), solver->kTot());

    c = 0;
    for (uint i = 0; i < solver->allPossibleReactions().size(); ++i)
    {
        CHECK_EQUAL(false, solver->isEmptyAddress(i));
        CHECK_EQUAL(Reaction::linearRateScale(), solver->allPossibleReactions().at(i)->rate());

        c += solver->allPossibleReactions().at(i)->rate();
        CHECK_EQUAL(c, solver->accuAllRates().at(i));
    }

    //activating a particle next to center particle. Should give 2 blocked reactions total.
    Site *neighbor = getBoxCenter(1);
    forceSpawnCenter(1);

    //spawning the neighbor should affect them both
    CHECK_EQUAL(2, SoluteParticle::affectedParticles().size());

    SoluteParticle::updateAffectedParticles();

    //However, the reaction blocked by center hasnt got a rate and is blocked, so it does not become available since it was never present.
    //Total extra possible reaction slots due to neighbor: 25
    //That is not the case for the reaction blocked by neighbor. So a total cavancy of one is expected had it not been for the fact that
    //another of the neighbors reactions fits in this slot.
    //Total extra possible reaction slots due to neighbor: 24
    //We expect a total vacancy of 0 and 24 + 26 stored reactions.
    CHECK_EQUAL(0, solver->availableReactionSlots().size());
    CHECK_EQUAL(50, solver->allPossibleReactions().size());
    CHECK_EQUAL(50, solver->accuAllRates().size());


    c = 0;

    for (uint i = 0; i < solver->allPossibleReactions().size(); ++i)
    {
        CHECK_EQUAL(false, solver->isEmptyAddress(i));

        c += solver->allPossibleReactions().at(i)->rate();
        CHECK_CLOSE(c, solver->accuAllRates().at(i), 0.00001);
    }


    //creating a with 26 fresh reactions.
    Site * distantCousin = getBoxCenter(0, 0, Site::nNeighborsLimit() + 1);

    SoluteParticle *p = new SoluteParticle();
    CHECK_EQUAL(true, solver->spawnParticle(p, NX()/2, NY()/2, NZ()/2 + Site::nNeighborsLimit() + 1, true));

    CHECK_EQUAL(1, SoluteParticle::affectedParticles().size());

    SoluteParticle::updateAffectedParticles();

    CHECK_EQUAL(0, solver->availableReactionSlots().size());
    CHECK_EQUAL(76, solver->allPossibleReactions().size());
    CHECK_EQUAL(76, solver->accuAllRates().size());

    //Now we move it so it kisses the center. This should give 2 vacant reactions.
    distantCousin->associatedParticle()->changePosition(NX()/2 - 1, NY()/2, NZ()/2);

    //dependent on the reach, the original neighbor is affected
    if (Site::nNeighborsLimit() > 1)
    {
        CHECK_EQUAL(3, SoluteParticle::affectedParticles().size());
    }
    else
    {
        CHECK_EQUAL(2, SoluteParticle::affectedParticles().size());
    }

    SoluteParticle::updateAffectedParticles();

    CHECK_EQUAL(2, solver->availableReactionSlots().size());
    CHECK_EQUAL(76, solver->allPossibleReactions().size());
    CHECK_EQUAL(76, solver->accuAllRates().size());

    //activating a new particle independent of the others should now only induce 24 more spots,
    //since two are already vacant.
    Site *distantCousin2 = getBoxCenter(0, 0, -(int)Site::nNeighborsLimit() - 1);

    SoluteParticle *p2 = new SoluteParticle();
    CHECK_EQUAL(true, solver->spawnParticle(p2, NX()/2, NY()/2, NZ()/2 - (int)Site::nNeighborsLimit() - 1, true));

    SoluteParticle::updateAffectedParticles();

    CHECK_EQUAL(0,  solver->availableReactionSlots().size());
    CHECK_EQUAL(100, solver->allPossibleReactions().size());
    CHECK_EQUAL(100, solver->accuAllRates().size());

    c = 0;
    for (uint i = 0; i < solver->allPossibleReactions().size(); ++i)
    {
        c += solver->allPossibleReactions().at(i)->rate();
        CHECK_CLOSE(c, solver->accuAllRates().at(i), 0.00001);
    }


    //deactivating the neighbor and the cousin should create two decoupled systems
    solver->despawnParticle(neighbor->associatedParticle());
    solver->despawnParticle(getBoxCenter(-1)->associatedParticle());

    SoluteParticle::updateAffectedParticles();

    //deactivating the two neighbors, which has 1 blocked reaction each, frees up 25 slots each. However, two of the original
    //center reactions now becomes active, and filles two slots, totalling 48 free slots after deactivation.
    CHECK_EQUAL(48,  solver->availableReactionSlots().size());
    CHECK_EQUAL(100, solver->allPossibleReactions().size());
    CHECK_EQUAL(100, solver->accuAllRates().size());

    CHECK_CLOSE(2*26*Reaction::linearRateScale(), solver->kTot(), 0.00001);

    c = 0;
    uint i = 0;
    for (SoluteParticle *particle : solver->particles())
    {
        particle->forEachActiveReactionDo([&] (Reaction *r)
        {
            if (solver->isEmptyAddress(i))
            {
                i++;
                return;
            }
            c += r->rate();
            CHECK_CLOSE(c, solver->accuAllRates().at(i), 0.00001);
            i++;

        });
    }

    //removing both the particles should even out everything to the maximum amount of reactions ever existing at once.

    solver->despawnParticle(center->associatedParticle());
    solver->despawnParticle(distantCousin2->associatedParticle());

    SoluteParticle::updateAffectedParticles();

    CHECK_EQUAL(100, solver->availableReactionSlots().size());
    CHECK_EQUAL(100, solver->allPossibleReactions().size());
    CHECK_EQUAL(100, solver->accuAllRates().size());

    CHECK_CLOSE(0, solver->kTot(), 0.00001);

    for (uint i = 0; i < solver->allPossibleReactions().size(); ++i)
    {
        CHECK_CLOSE(0, solver->accuAllRates().at(i), 0.00001);
    }


    Reaction::setLinearRateScale(1);

}

void testBed::testReactionShuffler()
{
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

    for (DummyReaction *r : allReacs)
    {
        delete r;
    }


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

    forceSpawnCenter();

    CHECK_EQUAL(ParticleStates::solvant, center->associatedParticle()->particleState());

    solver->despawnParticle(center->associatedParticle());

    for (int i = -1; i <= 1; ++i)
    {
        for (int j = -1; j <= 1; ++j)
        {
            for (int k = -1; k <= 1; ++k)
            {
                forceSpawnCenter(i, j, k);
            }
        }
    }

    forceSpawnCenter(3, 0, 0);

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
    solver->clearSites();

    Site::finalizeBoundaries();

    solver->setBoxSize({15, 15, 15});

    Site::resetNNeighborsLimitTo(5);

    solver->initializeSites();

    Site::initializeBoundaries();


    forceSpawnCenter();

    auto NNSUMBF = [&] (Site * site) -> uint
    {
        uint nn = 0;

        site->associatedParticle()->forEachNeighborSiteDo([&nn] (Site *neighbor)
        {
            if (neighbor->isActive())
            {
                nn++;
            }
        });

        return nn;
    };

    CHECK_EQUAL(0, NNSUMBF(getBoxCenter()));

    uint c = 0;
    uvec nn(Site::nNeighborsLimit(), fill::zeros);

    Site::forEachNeighborDo_sendPath(NX()/2, NY()/2, NZ()/2, [&] (Site *site, int dx, int dy, int dz)
    {
        (void) site;

        uint level = Site::levelMatrix(dx + Site::nNeighborsLimit(), dy + Site::nNeighborsLimit(), dz + Site::nNeighborsLimit());

        solver->forceSpawnParticle(Site::boundaries(0, 0)->transformCoordinate((int)NX()/2 + dx),
                                   Site::boundaries(1, 0)->transformCoordinate((int)NY()/2 + dy),
                                   Site::boundaries(2, 0)->transformCoordinate((int)NZ()/2 + dz));

        c++;
        nn(level)++;

        CHECK_EQUAL(c, NNSUMBF(getBoxCenter()));

    });

    Site::forEachNeighborDo_sendIndices(NX()/2, NY()/2, NZ()/2, [&] (Site *site, uint i, uint j, uint k)
    {
        uint level = Site::levelMatrix(i, j, k);

        solver->despawnParticle(site->associatedParticle());

        c--;
        nn(level)--;


        CHECK_EQUAL(c, NNSUMBF(getBoxCenter()));

    });

}

void testBed::testAccuAllRates()
{

    solver->initializeCrystal(0.3);

    solver->getRateVariables();

    uint NC = 100;
    Reaction *r;

    for (uint i = 0; i < NC; ++i)
    {
        r = solver->allPossibleReactions().at(KMC_RNG_UNIFORM()*solver->allPossibleReactions().size());
        r->execute();

        solver->getRateVariables();
    }

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

    activateAllSites();

    CHECK_EQUAL(NX()*NY()*NZ(), SoluteParticle::nParticles());

    solver->forEachSiteDo([&] (uint x, uint y, uint z, Site *_site)
    {

        site = _site;

        CHECK_EQUAL(x, site->associatedParticle()->x());
        CHECK_EQUAL(y, site->associatedParticle()->y());
        CHECK_EQUAL(z, site->associatedParticle()->z());

        CHECK_EQUAL(site, solver->getSite(site->associatedParticle()->x(),
                                          site->associatedParticle()->y(),
                                          site->associatedParticle()->z()));

    });


    for (uint i = 0; i < NX(); ++i)
    {

        for (uint j = 0; j < NX(); ++j)
        {

            for (uint k = 0; k < NX(); ++k)
            {
                site = solver->getSite(i, j, k);

                CHECK_EQUAL(i, site->associatedParticle()->x());
                CHECK_EQUAL(j, site->associatedParticle()->y());
                CHECK_EQUAL(k, site->associatedParticle()->z());
            }
        }
    }

    deactivateAllSites();

    forceNewBoundaries(Boundary::Periodic);

    solver->forEachSiteDo([&] (uint x, uint y, uint z, Site *_site)
    {
        site = _site; //hack to fix autocompletion.

        Boundary::setupCurrentBoundaries(x, y, z);

        for (const int &dx : Site::originTransformVector())
        {
            nx = (int)x + dx;
            xt = Boundary::currentBoundaries(0)->transformCoordinate(nx);

            for (const int &dy : Site::originTransformVector())
            {
                ny = (int)y + dy;
                yt = Boundary::currentBoundaries(1)->transformCoordinate(ny);

                for (const int &dz : Site::originTransformVector())
                {
                    nz = (int)z + dz;
                    zt = Boundary::currentBoundaries(2)->transformCoordinate(nz);

                    CHECK_EQUAL(false, Site::neighborhood(x, y, z, dx, dy, dz) == NULL);
                    CHECK_EQUAL(Site::neighborhood(x, y, z, dx, dy, dz), solver->getSite(nx, ny, nz));

                    if ((nx < 0 || nx >= (int)NX()) ||
                            (ny < 0 || ny >= (int)NY()) ||
                            (nz < 0 || nz >= (int)NZ()))
                    {
                        CHECK_EQUAL(*solver->getSite(xt, yt, zt), *solver->getSite(nx, ny, nz));
                    }

                }
            }
        }

    });

    forceNewBoundaries(Boundary::Edge);

    solver->forEachSiteDo([&] (uint x, uint y, uint z, Site *_site)
    {
        site = _site; //hack to fix autocompletion.

        Boundary::setupCurrentBoundaries(x, y, z);

        for (const int &dx : Site::originTransformVector())
        {
            nx = (int)x + dx;
            xt = Boundary::currentBoundaries(0)->transformCoordinate(nx);

            for (const int &dy : Site::originTransformVector())
            {
                ny = (int)y + dy;
                yt = Boundary::currentBoundaries(1)->transformCoordinate(ny);

                for (const int &dz : Site::originTransformVector())
                {
                    nz = (int)z + dz;
                    zt = Boundary::currentBoundaries(2)->transformCoordinate(nz);

                    CHECK_EQUAL(Site::neighborhood(x, y, z, dx, dy, dz), solver->getSite(nx, ny, nz));

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
