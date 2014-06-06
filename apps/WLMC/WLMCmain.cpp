#include <kMC>
#include <libconfig_utils/libconfig_utils.h>

#include <vector>

using namespace libconfig;
using namespace kMC;

using std::vector;

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

        cout << Seed::initialSeed << endl;


    }

protected:

    void execute()
    {

        moveParticles();

        DOS = normalise(DOS);

        if (nTimesExecuted() % 1000 != 0 || nTimesExecuted() == 0)
        {
            return;
        }

        cout << "N = " << SoluteParticle::nParticles() << " f=" << f << "  " << f/fCrit << endl;

        double flatness = estimateFlatness();

        output(flatness);

        if (flatness >= hCrit)
        {
            f = fReducer(f);
            resetCounts();

            if (f < fCrit)
            {
                prepNextOccupancyLevel();
            }
        }
    }

private:

    const uint nbins;

    vec maxEnergy;
    vec minEnergy;
    double energySpan;

    double meanFlatness;
    uint flatnessCounter;

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
    const uint nStart = 14;


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

        bool success = join_rows(join_rows(E, DOS), vc).eval().save(file.str());
        if (!success)
        {
            cout << "failed to dump data to file " << endl;
        }

    }

    double estimateFlatness()
    {
        uvec visitCounts_test = visitCounts;

        uint min = arma::min(visitCounts_test);

        double mean = 0;
        uint c = 0;
        for (const uint &vc : visitCounts_test)
        {
            if (vc != unsetCount)
            {
                mean += vc;
                c++;
            }
        }

        mean /= c;

        double flatness = min/mean;

        //        uint n = 100000;

        //        return (nTimesExecuted()%n)/double(n-1)*0.85;

        meanFlatness += flatness;
        flatnessCounter++;

        return meanFlatness/flatnessCounter;
    }

    void resetCounts()
    {
        visitCounts.fill(unsetCount);
        meanFlatness = 0;
        flatnessCounter = 0;
    }

    void initializeNewCycle()
    {
        resetCounts();

        DOS.ones();
        f = f0;

        energySpan = maxEnergy(count) - minEnergy(count);

        solver()->clearParticles();

        for (uint c = 0; c < count; ++c)
        {
            solver()->insertRandomParticle(0, false, false);
        }

        solver()->dumpLAMMPS(1337);

        KMCDebugger_Assert(energySpan, !=, 0);
    }

    double nMovesDist() const
    {
//        return 1 + SoluteParticle::nParticles()*KMC_RNG_UNIFORM();
//        return SoluteParticle::nParticles() - (SoluteParticle::nParticles()*(1 - std::sqrt(1-KMC_RNG_UNIFORM())) + 1);
//         return SoluteParticle::nParticles()*(1 - std::sqrt(1-KMC_RNG_UNIFORM())) + 1;
//        return SoluteParticle::nParticles()/2;
        return 1;
    }

    template<typename T>
    bool is_in(const vector<T> &_vec, const T &_val) const
    {
        return std::find(_vec.begin(), _vec.end(), _val) != _vec.end();
    }


    void moveParticles()
    {
        SoluteParticle* particle;

        uint which, where;

        uint xd = 0;
        uint yd = 0;
        uint zd = 0;

        uint nMoves = nMovesDist();

        vector<SoluteParticle*> selectedParticles;
        umat origPos(nMoves, 3);

        double E0 = SoluteParticle::totalEnergy();

        for (uint n = 0; n < nMoves; ++n)
        {
            which = KMC_RNG_UNIFORM()*(SoluteParticle::nParticles() - n);

            for (uint p = 0; p < which; ++p)
            {
                if (is_in(selectedParticles, solver()->particle(p)))
                {
                    which++;
                }
            }

            particle = solver()->particle(which);

            selectedParticles.push_back(particle);

            origPos(n, 0) = particle->x();
            origPos(n, 1) = particle->y();
            origPos(n, 2) = particle->z();


            where = (solver()->volume() - SoluteParticle::nParticles())*KMC_RNG_UNIFORM();

            findAvailableSite(where, xd, yd, zd);

            particle->changePosition(xd, yd, zd);

        }

        double E1 = SoluteParticle::totalEnergy();

        uint prevBin = getBin(E0);
        uint newBin = getBin(E1);

        double prevDOS = DOS(prevBin);
        double newDOS = DOS(newBin);

        bool accepted = true;

        if (prevDOS < newDOS)
        {
            accepted = KMC_RNG_UNIFORM() < prevDOS/newDOS;
        }

        if (accepted)
        {
            registerVisit(newBin);
        }

        else
        {
            registerVisit(prevBin);

            for (uint n = 0; n < nMoves; ++n)
            {
                uint nr = nMoves - n - 1;

                xd = origPos(nr, 0);
                yd = origPos(nr, 1);
                zd = origPos(nr, 2);

                selectedParticles.at(nr)->changePosition(xd, yd, zd);

            }

            KMCDebugger_AssertClose(E0, SoluteParticle::totalEnergy(), 1E-3, "reset failed.", SoluteParticle::totalEnergy());

        }

    }


    void findAvailableSite(const uint where, uint &xd, uint &yd, uint &zd)
    {
        uint search = 0;

        Site* site;

        for (uint x = 0; x < NX(); ++x)
        {
            for (uint y = 0; y < NY(); ++y)
            {
                for (uint z = 0; z < NZ(); ++z)
                {
                    site = solver()->getSite(x, y, z);

                    if (!site->isActive())
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
        solver()->forceSpawnParticle(NX()/2 + 1, NY()/2, NZ()/2);

        n = 2;
        do
        {
            parseForExtrema(+1);
            maxEnergy(n) = SoluteParticle::totalEnergy();

            n++;

            solver()->insertRandomParticle(0, false, false);

        } while (n < maxEnergy.n_elem);

    }



};



void initializeWLMC(KMCSolver *solver, const Setting &root)
{

    const Setting &initCFG = getSetting(root, "Initialization");

    const uint &nbins = getSetting<uint>(initCFG, "nbins");

    const double &f0 = getSetting<double>(initCFG, "f0");
    const double &fCrit = 1 + getSetting<double>(initCFG, "fCrit");
    const double &hCrit = getSetting<double>(initCFG, "hCrit");


    KMCEvent *wlmc = new WLMCEvent(nbins, f0, fCrit, hCrit);
    wlmc->setManualPriority(0);

    solver->addEvent(wlmc);

    uint solverEventAddress = solver->solverEvent()->getAddress();

    solver->mainLattice()->removeEvent(solverEventAddress);

}
