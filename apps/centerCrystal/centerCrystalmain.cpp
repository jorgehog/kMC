#include <kMC>
#include <libconfig_utils/libconfig_utils.h>

using namespace libconfig;
using namespace kMC;


void initialize_centerCrystal(KMCSolver * solver, const Setting & root);

int main()
{

    Config cfg;
    wall_clock t;


    cfg.readFile("infiles/centerCrystal.cfg");

    const Setting & root = cfg.getRoot();


    KMCDebugger_SetFilename("centerCrystal");

    KMCDebugger_SetEnabledTo(getSetting<int>(root, "buildTrace") == 0 ? false : true);

    KMCSolver::enableDumpLAMMPS(true);
    KMCSolver* solver = new KMCSolver(root);

    initialize_centerCrystal(solver, root);


    t.tic();

    solver->mainloop();

    cout << "Simulation ended after " << t.toc() << " seconds" << endl;


    KMCDebugger_DumpFullTrace();

    return 0;

}

class Sphericity : public KMCEvent
{
public:
    Sphericity() : KMCEvent("Sphericity", "", true, true) {}

protected:

    void execute()
    {
        uint A = 0;
        solver()->forEachSiteDo([&A] (uint x, uint y, uint z, Site* _site)
        {
            Site *site = _site;

            if (!site->isActive())
            {
                return;
            }

            if (site->associatedParticle()->isCrystal())
            {
                return;
            }

            if (!site->hasNeighboring(x, y, z, ParticleStates::crystal))
            {
               return;
            }

            uint nC = site->countNeighboring(x, y, z, ParticleStates::crystal);

            if (nC == 1)
            {
                A += 3;
            }

            else if (nC == 2 || nC == 3)
            {
                A += 2;
            }

            else
            {
                A += 1;
            }

        });

        setValue(pi3root*pow(6.0*(SoluteParticle::nCrystals() + A), 2./3)/A);
    }

private:

    static const double pi3root;

};

const double Sphericity::pi3root = pow(datum::pi, 1./3);


class TotalEnergy : public KMCEvent
{
public:

    TotalEnergy() : KMCEvent("TotalEnergy", "E* N", true, true) {}

protected:

    void execute()
    {
        setValue(SoluteParticle::totalEnergy()/SoluteParticle::nParticles());
    }

};

void initialize_centerCrystal(KMCSolver * solver, const Setting & root)
{

    solver->initializeCrystal(getSetting<double>(root, {"Initialization", "RelativeSeedSize"}));
    solver->initializeSolutionBath();
//    solver->initializeFromXYZ("/home/jorgen/code/build-kMC-Desktop_Qt_5_2_1_GCC_64bit-Release/apps/centerCrystal/outfiles", 37229);

    solver->addEvent(new Sphericity());
    solver->addEvent(new TotalEnergy());
}
