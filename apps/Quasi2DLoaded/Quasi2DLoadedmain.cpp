#include "quasidiffusionevents.h"
#include "quasidiffusion.h"

#include <commonkmcevents.h>


#include <libconfig_utils/libconfig_utils.h>

#include <HDF5Wrapper/include/hdf5wrapper.h>

using namespace libconfig;
using namespace kMC;


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

//    uint HMM = 10;
//    for (uint i= 0; i < solver->NX()/HMM; ++i)
//    {
//        heightmap(i*HMM) = KMC_RNG_UNIFORM()*1000;

//        for (uint j = 1; j < HMM; ++j)
//        {
//            heightmap(i*HMM + j) = heightmap(i*HMM);
//        }
//    }

    DumpHeighmap dumpHeightmap(heightmap);
    HeightRMS heightRMS(heightmap);
    AutoCorrHeight autocorr(heightmap);


    heightRMS.setDependency(dumpHeightmap);
    autocorr.setDependency(dumpHeightmap);

    solver->addEvent(new TotalTime());
    solver->addEvent(dumpHeightmap);
    solver->addEvent(heightRMS);
#ifndef NDEBUG
    cout << "checking rates." << endl;
    solver->addEvent(new RateChecker());
#endif

    const uint &h0 = getSetting<uint>(initCFG, "h0");
    const uint &therm = getSetting<uint>(initCFG, "therm");

    const double &Eb = getSetting<double>(initCFG, "Eb");
    const double &EsMax = getSetting<double>(initCFG, "EsMax")*Eb;
    const double &EsInit = getSetting<double>(initCFG, "EsInit")*Eb;

    const bool acf = getSetting<uint>(initCFG, "acf") == 1;

    const bool useWall = getSetting<uint>(initCFG, "useWall") == 1;

    const bool concEquil = getSetting<uint>(initCFG, "concEquil") == 1;
    const bool reset = getSetting<uint>(initCFG, "reset") == 1;
    const double &gCrit = getSetting<double>(initCFG, "gCrit");
    const double &treshold = getSetting<double>(initCFG, "treshold");
    const uint &N = getSetting<uint>(initCFG, "N");

    const uint &wallOnsetCycle = getSetting<uint>(initCFG, "wallOnsetCycle");

    solver->enableLocalUpdating(false);

    //Override standard diffusion. Necessary for quasi diffusive simulations.
    solver->setDiffusionType(KMCSolver::DiffusionTypes::None);

    MovingWall wallEvent(h0, EsMax, EsInit, heightmap);

    QuasiDiffusionSystem system(heightmap, Eb, wallEvent);

    if (acf)
    {
        solver->addEvent(autocorr);
    }

    NNeighbors nNeighbors(system);
    solver->addEvent(nNeighbors);

    Cumulant cumulant(system);
    solver->addEvent(cumulant);

    cumulant.setOnsetTime(therm);
    heightRMS.setOnsetTime(therm);
    autocorr.setOnsetTime(therm);

    EqConc eqC;
    ConcEquilibriator cc(eqC, N, gCrit, treshold);

    eqC.setDependency(nNeighbors);

    if (useWall)
    {
        eqC.setDependency(wallEvent);
        wallEvent.setDependency(solver->solverEvent());
        wallEvent.setOnsetTime(wallOnsetCycle);

        solver->addEvent(wallEvent);

        eqC.setOnsetTime(wallOnsetCycle + therm);
        cc.setOnsetTime(wallOnsetCycle + therm);
    }
    else
    {
        eqC.setOnsetTime(therm);
        cc.setOnsetTime(therm);
    }

    if (concEquil)
    {
        solver->addEvent(eqC);
        solver->addEvent(cc);
    }

    for (uint site = 0; site < solver->NX(); ++site)
    {
        SoluteParticle*  particle = solver->forceSpawnParticle(site, 0, 0);
        particle->addReaction(new LeftHopPressurized(particle, system));
        particle->addReaction(new RightHopPressurized(particle, system));
        particle->addReaction(new DepositionMirrorImageArhenius(particle, system));
        particle->addReaction(new Dissolution(particle, system));
    }

    H5Wrapper::Member &sizeMember = h5root.addMember(solver->NX());

    stringstream s;
    s << dynamic_cast<QuasiDiffusionReaction*>(solver->particle(0)->reactions().at(0))->numericDescription();
    s << "_n_" << runID;

    H5Wrapper::Member &potentialMember = sizeMember.addMember(s.str());

    if (concEquil && reset)
    {
        solver->mainloop();
        potentialMember.addData("eqConc", eqC.eqConc());

        solver->mainLattice()->removeEvent(&eqC);
        solver->mainLattice()->removeEvent(&cc);
    }

    t.tic();
    solver->mainloop();
    cout << "Simulation ended after " << t.toc() << " seconds" << endl;

    potentialMember.addData("heightmap", heightmap, overwrite);
    potentialMember.addData("ignisData", solver->mainLattice()->storedEventValues(), overwrite);
    potentialMember.addData("ignisEventDescriptions", solver->mainLattice()->outputEventDescriptions(), overwrite);
    potentialMember.addData("AutoCorr", autocorr.acf(), overwrite);

    if (concEquil && !reset)
    {
        potentialMember.addData("eqConc", eqC.eqConc());
    }

    return 0;

}

