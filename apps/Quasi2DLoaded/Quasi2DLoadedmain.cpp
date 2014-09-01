#include "quasidiffusion.h"
#include "quasidiffusionevents.h"

#include <commonkmcevents.h>


#include <libconfig_utils/libconfig_utils.h>

#include <HDF5Wrapper/include/hdf5wrapper.h>

using namespace libconfig;


ivec *initializeQuasi2DLoaded(KMCSolver * solver, const Setting & root, const uint l, const double beta);

int main()
{

    Config cfg;
    wall_clock t;


    cfg.readFile("infiles/Quasi2DLoaded.cfg");

    const Setting & root = cfg.getRoot();


    KMCDebugger_SetFilename("Quasi2DLoaded");

    KMCDebugger_SetEnabledTo(getSetting<int>(root, "buildTrace") == 0 ? false : true);

    KMCSolver::enableLocalUpdating(false);
    KMCSolver::enableDumpLAMMPS(false);
    MainLattice::enableEventFile(true, "ignisQuasi2Dloaded.arma", 1, 1);

    KMCSolver* solver = new KMCSolver(root);

    H5Wrapper::Root h5root("Quasi2D.h5");

    const Setting &initCFG = getSetting(root, "Initialization");

    const uint &lStart = getSetting<uint>(initCFG, "lStart");
    const uint &lEnd   = getSetting<uint>(initCFG, "lEnd");
    const uint &lStep  = getSetting<uint>(initCFG, "lStep");

    const double &tStart = getSetting<double>(initCFG, "tStart");
    const double &tEnd   = getSetting<double>(initCFG, "tEnd");
    const uint &nTemps  = getSetting<uint>(initCFG, "nTemps");

    vec temps = linspace(tStart, tEnd, nTemps);
    uint N = temps.size();

    t.tic();

    for (uint l = lStart; l <= lEnd; l += lStep)
    {
        for (uint i = 0; i < N; ++i)
        {
            double t = temps(i);

            cout << "Running l b = " << l << " " << t << endl;

            ivec* heightmap = initializeQuasi2DLoaded(solver, root, l, t);

            H5Wrapper::Member &sizeMember = h5root.addMember(l);

            H5Wrapper::Member &potentialMember = sizeMember.addMember(QuasiDiffusionReaction::potentialString());

            solver->mainloop();
            potentialMember.addData("heightmap", *heightmap);
            potentialMember.addData("ignisData", KMCEvent::eventMatrix());
            potentialMember.addData("ignisEventDescriptions", KMCEvent::outputEventDescriptions());

            potentialMember.file()->flush(H5F_SCOPE_GLOBAL);
            solver->reset();

            delete heightmap;
        }
    }

    cout << "Simulation ended after " << t.toc() << " seconds" << endl;

    KMCDebugger_DumpFullTrace();

    return 0;

}

ivec* initializeQuasi2DLoaded(KMCSolver *solver, const Setting &root, const uint l, const double beta)
{
    (void) root;

    solver->resetBoxSize(l, 1, 1);
    DiffusionReaction::setBeta(beta);

    ivec* heighmap = new ivec(l, fill::zeros);

    //Override standard diffusion. Necessary for quasi diffusive simulations.
    solver->setDiffusionType(KMCSolver::DiffusionTypes::None);

    //BAD PRATICE WITH POINTERS.. WILL FIX..
    MovingWall *wallEvent = new MovingWall(20, 10, 3, *heighmap);

    for (uint site = 0; site < l; ++site)
    {

        SoluteParticle* particle = solver->forceSpawnParticle(site, 0, 0);
        particle->addReaction(new LeftHop(particle, *heighmap, *wallEvent));
        particle->addReaction(new RightHop(particle, *heighmap, *wallEvent));
        particle->addReaction(new Deposition(particle, *heighmap, *wallEvent));
        particle->addReaction(new Dissolution(particle, *heighmap, *wallEvent));

    }

    QuasiDiffusionReaction::initialize();

    solver->addEvent(wallEvent);
    solver->addEvent(new DumpHeighmap(*heighmap));
    solver->addEvent(new TotalTime());
    solver->addEvent(new heightRMS(*heighmap));

    return heighmap;

}
