#include <kMC>
#include <libconfig_utils/libconfig_utils.h>

using namespace libconfig;
using namespace kMC;


void initializeWLMC(KMCSolver * solver, const Setting & root);

int main()
{

    Config cfg;
    wall_clock t;


    cfg.readFile("infiles/WLMC.cfg");

    const Setting & root = cfg.getRoot();


    KMCDebugger_SetFilename("WLMC");

    KMCDebugger_SetEnabledTo(getSetting<int>(root, "buildTrace") == 0 ? false : true);

    KMCSolver::enableDumpLAMMPS(false);
    KMCSolver::enableDumpXYZ(false);
    MainLattice::enableEventFile(false);

    KMCSolver* solver = new KMCSolver(root);

    initializeWLMC(solver, root);


    t.tic();

    solver->mainloop();

    cout << "Simulation ended after " << t.toc() << " seconds" << endl;



    KMCDebugger_DumpFullTrace();

    return 0;

}


class WLMCEvent : public KMCEvent
{
public:

    WLMCEvent(const uint nbins) :
        KMCEvent("kMC::WLMC", "", true),
        nbins(nbins)
    {

    }

    void initialize()
    {

        minEnergy.set_size(solver()->volume());
        maxEnergy.set_size(solver()->volume());

        findAllExtrema();

        count = nStart;

        visitCounts.set_size(nbins);
        DOS.set_size(nbins);

        initializeNewCycle();


    }

protected:

    void execute()
    {

        moveParticle();

        cout << "N=" << SoluteParticle::nParticles() << " " << count << endl;
        double flatness = estimateFlatness();

        if (flatness > 0.85)
        {
            f = sqrt(f);
            visitCounts.zeros();

            if (f < 1.000001)
            {
                prepNextOccupancyLevel();
            }
        }

        setValue(flatness);

        stringstream file;
        file << solver()->filePath() << "/stateDensity" << SoluteParticle::nParticles() << ".arma";
        vec E = linspace<vec>(minEnergy(count), maxEnergy(count), nbins);
        join_rows(join_rows(E(span(nbins/6, (nbins*5)/6)), DOS(span(nbins/6, (nbins*5)/6))), visitCounts(span(nbins/6, (nbins*5)/6))).eval().save(file.str());


    }

private:

    const uint nbins;

    vec maxEnergy;
    vec minEnergy;
    double energySpan;

    double f;

    uint count;

    vec visitCounts;
    vec DOS;

    const uint nSkipped = 5;
    const uint nStart = 10;


    double estimateFlatness()
    {
        vec visitSpan = visitCounts(span(nbins/3, (nbins*2)/3));

        double M = max(abs(visitSpan - mean(visitSpan)));

        double std = stddev(visitSpan);

        double flatness = std/M;

        cout << std << " " << M << " " << flatness << endl;

        return flatness;
    }

    void initializeNewCycle()
    {
        visitCounts.zeros();

        DOS.ones();
        f = datum::e;

        energySpan = maxEnergy(count) - minEnergy(count);

        KMCDebugger_Assert(energySpan, !=, 0);
    }

    void moveParticle()
    {
        uint which = KMC_RNG_UNIFORM()*SoluteParticle::nParticles();
        uint where = KMC_RNG_UNIFORM()*(solver()->volume() - SoluteParticle::nParticles() - 1);

        SoluteParticle* particle = solver()->particle(which);

        uint xd = 0;
        uint yd = 0;
        uint zd = 0;

        findAvailableSite(where, xd, yd, zd);

        performTrialMove(particle, xd, yd, zd);

    }

    void findAvailableSite(uint where, uint &xd, uint &yd, uint &zd)
    {
        uint search = 0;

        for (uint x = 0; x < NX(); ++x)
        {
            for (uint y = 0; y < NY(); ++y)
            {
                for (uint z = 0; z < NZ(); ++z)
                {
                    if (!solver()->getSite(x, y, z)->isActive())
                    {
                        if (search == where)
                        {
                            xd = x;
                            yd = y;
                            zd = z;

                            return;
                        }

                        search++;
                    }
                }
            }
        }

    }

    void performTrialMove(SoluteParticle *particle, uint xd, uint yd, uint zd)
    {

        double eNew = SoluteParticle::totalEnergy() + getEnergyDifference(particle, xd, yd, zd);

        uint prevBin = getBin(SoluteParticle::totalEnergy());
        uint newBin = getBin(eNew);

        double prevDOS = DOS(prevBin);
        double newDOS = DOS(newBin);

        if (KMC_RNG_UNIFORM() < prevDOS/newDOS)
        {
            particle->changePosition(xd, yd, zd);

            DOS(newBin) *= f;
            visitCounts(newBin)++;
        }

        else
        {
            DOS(prevBin) *= f;
            visitCounts(prevBin)++;
        }

    }

    uint getBin(double energy)
    {
        return (nbins - 1)*(energy - minEnergy(count))/energySpan;
    }

    void prepNextOccupancyLevel()
    {
        if (SoluteParticle::nParticles() == NX()*NY()*NZ() - nSkipped + 2)
        {
            solver()->exit();
        }

        bool spawned = false;
        solver()->forEachSiteDo([&spawned, this] (uint x, uint y, uint z, Site *site)
        {
            if (!site->isActive() && !spawned)
            {
                spawned = true;
                solver()->forceSpawnParticle(x, y, z);
            }
        });

        count++;

        initializeNewCycle();

    }


    void parseForExtrema(int sign)
    {

        double localExtrema;
        uint xd, yd, zd;

        bool done;

        do {

            done = true;

            if (sign == -1)
            {
                solver()->sortParticles([] (SoluteParticle *p1, SoluteParticle *p2) {return p1->energy() > p2->energy();});
            }
            else
            {
                solver()->sortParticles([] (SoluteParticle *p1, SoluteParticle *p2) {return p1->energy() < p2->energy();});
            }

            for (SoluteParticle *particle : solver()->particles())
            {
                localExtrema = 0;

                solver()->forEachSiteDo([&] (uint x, uint y, uint z, Site* site)
                {
                    if (site->isActive())
                    {
                        return;
                    }

                    double dE = getEnergyDifference(particle, x, y, z);

                    if (sign == 1 ? dE > localExtrema : dE < localExtrema)
                    {
                        localExtrema = dE;

                        xd = x;
                        yd = y;
                        zd = z;

                    }

                });

                if (localExtrema == 0)
                {
                    continue;
                }

                done = false;

                particle->changePosition(xd, yd, zd);

            }

        } while (!done);


        solver()->dumpLAMMPS(2*(SoluteParticle::nParticles() - 2) + (1 + sign)/2);

    }

    double getEnergyDifference(SoluteParticle *particle, const uint xd, const uint yd, const uint zd)
    {

        double eNew = 0;

        Site::forEachNeighborDo_sendIndices(xd, yd, zd, [&particle, &eNew] (Site *neighbor, uint i, uint j, uint k)
        {
            if (neighbor->isActive())
            {
                if (neighbor->associatedParticle() != particle)
                {
                    eNew += DiffusionReaction::potential(i, j, k);
                }
            }
        });

        return eNew - particle->energy();
    }


    void findAllExtrema()
    {
        uint n;

        solver()->forceSpawnParticle(0, 0, 0);
        solver()->forceSpawnParticle(NX() - 1, NY() - 1, NZ() - 1);

        n = 2;
        do
        {
            parseForExtrema(-1);
            minEnergy(n) = SoluteParticle::totalEnergy();

            n++;

            solver()->insertRandomParticle(0, false, false);

        } while (n < minEnergy.n_elem);

        solver()->clearParticles();



        solver()->forceSpawnParticle(NX()/2, NY()/2, NZ()/2);
        solver()->forceSpawnParticle(NX()/2, NY()/2, NZ()/2 + 1);

        n = 2;
        do
        {
            parseForExtrema(+1);
            maxEnergy(n) = SoluteParticle::totalEnergy();

            n++;

            solver()->insertRandomParticle(0, false, false);

        } while (n < maxEnergy.n_elem);

        solver()->clearParticles();

        for (uint c = 0; c < nStart; ++c)
        {
            solver()->insertRandomParticle();
        }
    }



};



void initializeWLMC(KMCSolver *solver, const Setting &root)
{

    const Setting &initCFG = getSetting(root, "Initialization");

    const uint &nbins = getSetting<uint>(initCFG, "nbins");


    KMCEvent *wlmc = new WLMCEvent(nbins);
    wlmc->setManualPriority(0);

    solver->addEvent(wlmc);

    uint solverEventAddress = solver->solverEvent()->getAddress();

    solver->mainLattice()->removeEvent(solverEventAddress);

}
