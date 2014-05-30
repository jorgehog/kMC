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

        minEnergy.set_size(solver()->volume() - 3);
        maxEnergy.set_size(solver()->volume() - 3);

        findAllExtrema();

        count = 0;

        visitCounts.set_size(nbins);
        DOS.set_size(nbins);

        resetHistograms();

        uint bin;
        double prevDOS, newDOS;
        updateHistograms(SoluteParticle::totalEnergy(), prevDOS, newDOS, bin);


    }

protected:

    void execute()
    {

        moveParticle();

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

    }

private:

    const uint nbins;
    const uint nCyclesPerPopulation = 100000;

    vec maxEnergy;
    vec minEnergy;
    double energySpan;

    double f;

    uint count;

    vec visitCounts;
    vec DOS;


    double estimateFlatness()
    {
        double m = mean(visitCounts);

        return max(abs(visitCounts - m))/max(visitCounts);
    }

    void resetHistograms()
    {
        visitCounts.zeros();
        DOS.ones();

        f = datum::e;

        energySpan = maxEnergy(count) - minEnergy(count);

        KMCDebugger_Assert(energySpan, !=, 0);
    }

    void updateHistograms(const double energy, double &prevDOS, double &newDOS, uint &bin)
    {
        bin = (nbins - 1)*(energy - minEnergy(count))/energySpan;

        visitCounts(bin)++;

        prevDOS = DOS(bin);

        DOS(bin) *= f;

        newDOS = DOS(bin);

    }

    void undoHistogram(const double prev, const uint bin)
    {
        DOS(bin) = prev;
        visitCounts(bin)--;
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

        double prevDOS, newDOS;
        uint bin;

        updateHistograms(eNew, prevDOS, newDOS, bin);


        if (KMC_RNG_UNIFORM() < prevDOS/newDOS)
        {
            particle->changePosition(xd, yd, zd);
        }

        else
        {
            undoHistogram(prevDOS, bin);

            updateHistograms(SoluteParticle::totalEnergy(), prevDOS, newDOS, bin);
        }

    }


    void prepNextOccupancyLevel()
    {
        stringstream file;
        file << solver()->filePath() << "/stateDensity" << SoluteParticle::nParticles() << ".arma";
        join_rows(DOS, visitCounts).eval().save(file.str());

        if (SoluteParticle::nParticles() == NX()*NY()*NZ() - 1)
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

        resetHistograms();

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
        uint c;

        solver()->forceSpawnParticle(0, 0, 0);
        solver()->forceSpawnParticle(NX() - 1, NY() - 1, NZ() - 1);

        c = 0;
        do
        {
            parseForExtrema(-1);
            minEnergy(c) = SoluteParticle::totalEnergy();

            c++;

            solver()->insertRandomParticle(0, false, false);

        } while (c < minEnergy.n_elem);

        solver()->clearParticles();



        solver()->forceSpawnParticle(NX()/2, NY()/2, NZ()/2);
        solver()->forceSpawnParticle(NX()/2, NY()/2, NZ()/2 + 1);

        c = 0;
        do
        {
            parseForExtrema(+1);
            maxEnergy(c) = SoluteParticle::totalEnergy();

            c++;

            solver()->insertRandomParticle(0, false, false);

        } while (c < maxEnergy.n_elem);

        solver()->clearParticles();


        solver()->insertRandomParticle();
        solver()->insertRandomParticle();
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
