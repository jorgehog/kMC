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

    WLMCEvent(const uint nbins, const double f0 = datum::e, const double fCrit=1.000001, const double hCrit = 0.85) :
        KMCEvent("kMC::WLMC", "", true),
        nbins(nbins),
        f0(f0),
        fCrit(fCrit),
        hCrit(hCrit)
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

        DOS = normalise(DOS);

        double flatness = estimateFlatness();

        output(flatness);

        if (flatness >= hCrit)
        {
            f = fReducer(f);
            resetCounts();

            if (f < fCrit)
            {
                prepNextOccupancyLevel();
                exit(1);
            }
        }

        cout << min(visitCounts) << "  " << mean(visitCounts) << endl;

    }

private:

    const uint nbins;

    vec maxEnergy;
    vec minEnergy;
    double energySpan;

    double f;
    double f0;
    double fCrit;
    double hCrit;

    function<double(double)> fReducer = [] (double fPrev) {return sqrt(fPrev);};

    uint count;


    uvec visitCounts;
    vec DOS;

    static constexpr uint unsetCount = std::numeric_limits<uint>::max();

    const uint nSkipped = 1;
    const uint nStart = 5;


    void output(double flatness)
    {
        setValue(flatness);

        stringstream file;
        file << solver()->filePath() << "/stateDensity" << SoluteParticle::nParticles() << ".arma";
        vec E = linspace<vec>(minEnergy(count), maxEnergy(count), nbins);

        vec vc = conv_to<vec>::from(visitCounts);

        double mean = 0;
        uint c = 0;
        for (uint i = 0; i < vc.n_elem; ++i)
        {
            if (vc(i) != unsetCount)
            {
                mean += vc(i);
                c++;
            }
        }

        mean /= c;

        for (uint i = 0; i < vc.n_elem; ++i)
        {
            if (vc(i) == unsetCount)
            {
                vc(i) = mean;
            }
        }

        join_rows(join_rows(E, DOS), vc).eval().save(file.str());

    }

    double estimateFlatness()
    {
        uint MAX = max(visitCounts);
        uint MIN = min(visitCounts);

        double span = MAX - MIN;

        double m = mean(visitCounts);



        double flatness = 1 - span/m;

        cout << uvec(find(visitCounts == unsetCount)).t() << endl;

        if (flatness < 0)
        {
            return 0;
        }
        uint n = 100000;

        return (nTimesExecuted()%n)/double(n-1)*0.85;

        return flatness;
    }

    void resetCounts()
    {
        visitCounts.fill(unsetCount);
    }

    void initializeNewCycle()
    {
        resetCounts();

        DOS.ones();
        f = f0;

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

        cout << SoluteParticle::totalEnergy() << endl;
        cout << maxEnergy(count) << endl;

        double prevDOS = DOS(prevBin);
        double newDOS = DOS(newBin);

        bool accepted = true;

        if (prevDOS < newDOS)
        {
            accepted = KMC_RNG_UNIFORM() < prevDOS/newDOS;
        }

        if (accepted)
        {
            particle->changePosition(xd, yd, zd);

            registerVisit(newBin);

        }

        else
        {
            registerVisit(prevBin);
        }

    }

    void registerVisit(const uint bin)
    {
        if (visitCounts(bin) == unsetCount)
        {
            visitCounts(bin) = 1;
        }
        else
        {
            visitCounts(bin)++;
        }

        DOS(bin) *= f;
    }

    uint getBin(double energy)
    {
        if (fabs(energy - maxEnergy(count)) < 1E-10)
        {
            return nbins - 1;
        }

        return nbins*(energy - minEnergy(count))/energySpan;
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

    const double &f0 = getSetting<double>(initCFG, "f0");
    const double &fCrit = getSetting<double>(initCFG, "fCrit");
    const double &hCrit = getSetting<double>(initCFG, "hCrit");


    KMCEvent *wlmc = new WLMCEvent(nbins, f0, fCrit, hCrit);
    wlmc->setManualPriority(0);

    solver->addEvent(wlmc);

    uint solverEventAddress = solver->solverEvent()->getAddress();

    solver->mainLattice()->removeEvent(solverEventAddress);

}
