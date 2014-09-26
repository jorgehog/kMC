#include "quasidiffusionevents.h"
#include "quasidiffusion.h"

#include <commonkmcevents.h>


#include <libconfig_utils/libconfig_utils.h>

#include <HDF5Wrapper/include/hdf5wrapper.h>

using namespace libconfig;
using namespace kMC;

void initializeQuasi2DLoaded(KMCSolver * solver, const Setting & root, ivec *heigthmap);

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

    const bool overwrite = getSetting<uint>(initCFG, "overwrite") == 1;


    ivec heightmap(solver->NX(), fill::zeros);

    DumpHeighmap dumpHeightmap(heightmap);
    HeightRMS heightRMS(heightmap);
    AutoCorrHeight autocorr(heightmap);


    heightRMS.setDependency(dumpHeightmap);
    autocorr.setDependency(dumpHeightmap);

    solver->addEvent(new TotalTime());
    solver->addEvent(dumpHeightmap);
    solver->addEvent(heightRMS);
    solver->addEvent(autocorr);
#ifndef NDEBUG
    cout << "checking rates." << endl;
    solver->addEvent(new RateChecker());
#endif

    initializeQuasi2DLoaded(solver, initCFG, &heightmap);

    H5Wrapper::Member &sizeMember = h5root.addMember(solver->NX());

    stringstream s;
    s << dynamic_cast<QuasiDiffusionReaction*>(solver->particle(0)->reactions().at(0))->numericDescription();
    s << "_n_" << runID;

    H5Wrapper::Member &potentialMember = sizeMember.addMember(s.str());

    t.tic();
    solver->mainloop();
    cout << "Simulation ended after " << t.toc() << " seconds" << endl;

    potentialMember.addData("heightmap", heightmap, overwrite);
    potentialMember.addData("ignisData", solver->mainLattice()->storedEventValues(), overwrite);
    potentialMember.addData("ignisEventDescriptions", solver->mainLattice()->outputEventDescriptions(), overwrite);
    potentialMember.addData("AutoCorr", autocorr.acf(), overwrite);

    potentialMember.file()->flush(H5F_SCOPE_GLOBAL);

    return 0;

}


void initializeQuasi2DLoaded(KMCSolver *solver, const Setting &initCFG, ivec *heigthmap)
{

    const uint &h0 = getSetting<uint>(initCFG, "h0");

    const double &Eb = getSetting<double>(initCFG, "Eb");
    const double &EsMax = getSetting<double>(initCFG, "EsMax")*Eb;
    const double &EsInit = getSetting<double>(initCFG, "EsInit")*Eb;

    const bool useWall = getSetting<uint>(initCFG, "useWall") == 1;
    const uint &wallOnsetCycle = getSetting<uint>(initCFG, "wallOnsetCycle");

    solver->enableLocalUpdating(false);

    //Override standard diffusion. Necessary for quasi diffusive simulations.
    solver->setDiffusionType(KMCSolver::DiffusionTypes::None);

    //BAD PRATICE WITH POINTERS.. WILL FIX..
    MovingWall *wallEvent = new MovingWall(h0, EsMax, EsInit, *heigthmap);

    wallEvent->setDependency(solver->solverEvent());
    wallEvent->setOnsetTime(wallOnsetCycle);

    if (useWall)
    {
        solver->addEvent(wallEvent);
    }

    for (uint site = 0; site < solver->NX(); ++site)
    {
        SoluteParticle* particle = solver->forceSpawnParticle(site, 0, 0);
        particle->addReaction(new LeftHopPressurized(particle, *heigthmap, Eb, *wallEvent));
        particle->addReaction(new RightHopPressurized(particle, *heigthmap, Eb, *wallEvent));
        particle->addReaction(new DepositionMirrorImageArhenius(particle, *heigthmap, Eb, *wallEvent));
        particle->addReaction(new Dissolution(particle, *heigthmap, Eb, *wallEvent));
    }


}
