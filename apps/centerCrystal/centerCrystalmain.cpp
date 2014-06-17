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

//    const uint &NX = solver->NX();
//    const uint &NY = solver->NY();
//    const uint &NZ = solver->NZ();

//    uint x, y, z1, z2;

//    vec thetas = linspace(0, datum::pi/2, 500);
//    vec phis   = linspace(0, 2*datum::pi, 500);

//    for (const double &theta : thetas)
//    {
//        for (const double &phi : phis)
//        {
//            x  = (NX - 1.0)/2.0 + round(NX*sin(theta)*cos(phi)/2.0);
//            y  = (NY - 1.0)/2.0 + round(NY*sin(theta)*sin(phi)/2.0);
//            z1 = (NZ - 1.0)/2.0 + round(NZ*cos(theta)         /2.0);
//            z2 = (NZ - 1.0)/2.0 - round(NZ*cos(theta)         /2.0);

//            if (!solver->getSite(x, y, z1)->isActive())
//            {
//                solver->forceSpawnParticle(x, y, z1);
//            }

//            if (!solver->getSite(x, y, z2)->isActive())
//            {
//                solver->forceSpawnParticle(x, y, z2);
//            }

//        }
//    }

//    solver->dumpLAMMPS(1337);


}
