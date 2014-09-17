#include "quasidiffusion.h"
#include "quasidiffusionevents.h"

#include <commonkmcevents.h>


#include <libconfig_utils/libconfig_utils.h>

#include <HDF5Wrapper/include/hdf5wrapper.h>

using namespace libconfig;

ivec *initializeQuasi2DLoaded(KMCSolver * solver, const Setting & root);

int main()
{

    Config cfg;
    wall_clock t;


    cfg.readFile("infiles/Quasi2DLoaded.cfg");

    const Setting & root = cfg.getRoot();


    KMCDebugger_SetFilename("Quasi2DLoaded");

    KMCDebugger_SetEnabledTo(getSetting<int>(root, "buildTrace") == 0 ? false : true);

//    KMCSolver::enableLocalUpdating(false);
    KMCSolver::enableDumpLAMMPS(false);

    KMCSolver* solver = new KMCSolver(root);

    string ignisOutputName = "ignisQuasi2Dloaded.ign";
    solver->mainLattice()->enableEventValueStorage(true, true, ignisOutputName, solver->filePath(), 1);

    H5Wrapper::Root h5root("Quasi2D.h5");

    const Setting &initCFG = getSetting(root, "Initialization");

    const uint &runID = getSetting<uint>(initCFG, "runID");

        t.tic();

        ivec* heightmap = initializeQuasi2DLoaded(solver, initCFG);

        H5Wrapper::Member &sizeMember = h5root.addMember(solver->NX());


        stringstream s;
        s << dynamic_cast<QuasiDiffusionReaction*>(solver->particle(0)->reactions().at(0))->numericDescription();
        s << "_n_" << runID;

        H5Wrapper::Member &potentialMember = sizeMember.addMember(s.str());

        solver->mainloop();

        potentialMember.addData("heightmap", *heightmap);
        potentialMember.addData("ignisData", solver->mainLattice()->storedEventValues());
        potentialMember.addData("ignisEventDescriptions", solver->mainLattice()->outputEventDescriptions());

        potentialMember.file()->flush(H5F_SCOPE_GLOBAL);
        solver->reset();

        delete heightmap;

        cout << "Simulation ended after " << t.toc() << " seconds" << endl;


    KMCDebugger_DumpFullTrace();

    return 0;

}


ivec* initializeQuasi2DLoaded(KMCSolver *solver, const Setting &initCFG)
{

    const uint &h0 = getSetting<uint>(initCFG, "h0");

    const double &Eb = getSetting<double>(initCFG, "Eb");
    const double &EsMax = getSetting<double>(initCFG, "EsMax")*Eb;
    const double &EsInit = getSetting<double>(initCFG, "EsInit")*Eb;

    ivec* heighmap = new ivec(solver->NX(), fill::zeros);

    solver->enableLocalUpdating(false);

    //Override standard diffusion. Necessary for quasi diffusive simulations.
    solver->setDiffusionType(KMCSolver::DiffusionTypes::None);

    //BAD PRATICE WITH POINTERS.. WILL FIX..
    MovingWall *wallEvent = new MovingWall(h0, EsMax, EsInit, *heighmap);

    for (uint site = 0; site < solver->NX(); ++site)
    {
        SoluteParticle* particle = solver->forceSpawnParticle(site, 0, 0);
        particle->addReaction(new LeftHopPressurized(particle, *heighmap, Eb, *wallEvent));
        particle->addReaction(new RightHopPressurized(particle, *heighmap, Eb, *wallEvent));
        particle->addReaction(new DepositionMirrorImageArhenius(particle, *heighmap, Eb, *wallEvent));
        particle->addReaction(new Dissolution(particle, *heighmap, Eb, *wallEvent));

    }

    solver->addEvent(wallEvent);
    solver->addEvent(new DumpHeighmap(*heighmap));
    solver->addEvent(new TotalTime());
    solver->addEvent(new heightRMS(*heighmap));

    solver->addEvent(new RateChecker());

    return heighmap;

}
