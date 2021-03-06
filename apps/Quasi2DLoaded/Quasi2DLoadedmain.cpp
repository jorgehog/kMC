#include "quasidiffusionevents.h"
#include "quasidiffusion.h"

#include <commonkmcevents.h>


#include <libconfig_utils/libconfig_utils.h>

#include <HDF5Wrapper/include/hdf5wrapper.h>

using namespace libconfig;
using namespace kMC;

string addProcEnding(string filename, string ending, string procending)
{
    stringstream s;
    s << filename << procending << "." << ending;

    return s.str();
}

int main(int argv, char** argc)
{

    string tail;

    if (argv == 2)
    {
        int proc = atoi(argc[1]);

        stringstream ending;
        ending << "_" << proc;

        tail = ending.str();
    }
    else
    {
        tail = "";
    }

    Config cfg;
    wall_clock t;

    string cfgName = "infiles/" + addProcEnding("Quasi2DLoaded", "cfg", tail);
    cfg.readFile(cfgName.c_str());

    const Setting & root = cfg.getRoot();


    KMCDebugger_SetFilename("Quasi2DLoaded");

    KMCDebugger_SetEnabledTo(getSetting<int>(root, "buildTrace") == 0 ? false : true);


    KMCSolver::enableDumpLAMMPS(false);

    KMCSolver* solver = new KMCSolver(root);

    string ignisOutputName = addProcEnding("ignisQuasi2Dloaded", "ign", tail);
    solver->mainLattice()->enableEventValueStorage(true, true, ignisOutputName, solver->filePath(), 1000);

    H5Wrapper::Root h5root(solver->filePath() +  addProcEnding("Quasi2D", "h5", tail));

    const Setting &initCFG = getSetting(root, "Initialization");

    const uint &runID = getSetting<uint>(initCFG, "runID");

    const bool overwrite = getSetting<uint>(initCFG, "overwrite") == 1;


    ivec heightmap(solver->NX(), fill::zeros);

//    uint HMM = 1;
//    for (uint i= 0; i < solver->NX()/HMM; ++i)
//    {
//        heightmap(i*HMM) = KMC_RNG_UNIFORM()*10;

//        for (uint j = 1; j < HMM; ++j)
//        {
//            heightmap(i*HMM + j) = heightmap(i*HMM);
//        }
//    }



    const uint &therm = getSetting<uint>(initCFG, "therm");

    const double &alpha = getSetting<double>(initCFG, "alpha");
    const double &mu = getSetting<double>(initCFG, "mu");

    const double &concAdd = getSetting<double>(initCFG, "concAdd");


    const double &sigma0 = getSetting<double>(initCFG, "sigma0");
    const double &r0 = getSetting<double>(initCFG, "r0");
    const double &Ew0dL = getSetting<double>(initCFG, "Ew0dL");
    const double Ew0 = -Ew0dL*heightmap.size();

    const bool acf = getSetting<uint>(initCFG, "acf") == 1;

    const uint &concEquilInt = getSetting<uint>(initCFG, "concEquil");
    const uint &shadowingInt = getSetting<uint>(initCFG, "shadowing");
    const uint &useDiffusionInt = getSetting<uint>(initCFG, "useDiffusion");
    const uint &isotropicDiffusionInt = getSetting<uint>(initCFG, "isotropicDiffusion");
    const uint &useWallInt = getSetting<uint>(initCFG, "useWall");
    const uint &resetInt = getSetting<uint>(initCFG, "reset");

    const bool useWall = useWallInt == 1;
    const bool useDiffusion = useDiffusionInt == 1;
    const bool useIsotropicDiffusion = isotropicDiffusionInt == 1;
    const bool useShadowing = shadowingInt == 1;

    const bool concEquil = concEquilInt == 1;
    const bool reset = resetInt == 1;
    const uint &N = getSetting<uint>(initCFG, "N");
    const uint &nRounds = getSetting<uint>(initCFG, "nRounds");

    const uint &wallOnsetCycle = getSetting<uint>(initCFG, "wallOnsetCycle");

    solver->enableLocalUpdating(false);

    //Override standard diffusion. Necessary for quasi diffusive simulations.
    solver->setDiffusionType(KMCSolver::DiffusionTypes::None);


    MovingWall wallEvent(Ew0, sigma0, r0, heightmap);

    QuasiDiffusionSystem system(heightmap, wallEvent, alpha, mu);

    DumpHeighmap dumpHeightmap(heightmap);
    solver->addEvent(dumpHeightmap);

    NNeighbors nNeighbors(system);
    solver->addEvent(nNeighbors);

    EqConc eqC(useShadowing);
    ConcEquilibriator cc(system, eqC, nRounds, N);

    eqC.setDependency(nNeighbors);

    if (useWall)
    {
        wallEvent.setDependency(dumpHeightmap);
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

#ifndef NDEBUG
    cout << "checking rates." << endl;
    RateChecker rateChecker;
    solver->addEvent(rateChecker);
#endif

    for (uint site = 0; site < solver->NX(); ++site)
    {
        SoluteParticle*  particle = solver->forceSpawnParticle(site, 0, 0);
        if (useDiffusion)
        {
            if (useIsotropicDiffusion)
            {
                particle->addReaction(new LeftHopIsotropic(particle, system));
                particle->addReaction(new RightHopIsotropic(particle, system));
            }
            else
            {
                particle->addReaction(new LeftHopDownOnly(particle, system));
                particle->addReaction(new RightHopDownOnly(particle, system));
            }
        }

        if (!useShadowing)
        {
            particle->addReaction(new DepositionMirrorImageArheniusNoShadowing(particle, system));
        }
        else
        {
            particle->addReaction(new DepositionMirrorImageArhenius(particle, system));
        }

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

        if (!cc.finalized())
        {
            cc.finalizeAverages();
        }

        potentialMember.addData("muEq", cc.averageMu(), overwrite);
        potentialMember.addData("muEqError", cc.error(), overwrite);

        system.setMu(cc.averageMu());

        solver->mainLattice()->removeEvent(&eqC);
        solver->mainLattice()->removeEvent(&cc);

        for (SoluteParticle *particle : solver->particles())
        {
            for (Reaction *reaction : particle->reactions())
            {
                if (reaction->isAllowed())
                {
                    reaction->registerUpdateFlag(QuasiDiffusionReaction::UpdateFlags::CALCULATE);
                }
            }
        }

        sleep(3);
    }

    if (concAdd != 0)
    {
        system.setMu(concAdd + system.mu());
    }

    HeightRMS heightRMS(heightmap);
    AutoCorrHeight autocorr(heightmap);


    heightRMS.setDependency(dumpHeightmap);
    autocorr.setDependency(dumpHeightmap);

    solver->addEvent(new TotalTime());
    solver->addEvent(heightRMS);

    if (acf)
    {
        solver->addEvent(autocorr);
    }


    Cumulant cumulant(system);
//    solver->addEvent(cumulant);

    SurfaceSize size(system);
    solver->addEvent(size);

    SurfaceSizeLocal sizeLocal;
    sizeLocal.setDependency(size);
    solver->addEvent(sizeLocal);

    autocorr.setOnsetTime(therm);
    nNeighbors.setOnsetTime(therm);
    cumulant.setOnsetTime(therm);
    size.setOnsetTime(therm);
    sizeLocal.setOnsetTime(therm);

    t.tic();
    solver->mainloop();
    cout << "Simulation ended after " << t.toc() << " seconds" << endl;


    potentialMember.addData("usewall", useWallInt, overwrite);
    potentialMember.addData("sigma0", sigma0, overwrite);
    potentialMember.addData("r0", r0, overwrite);
    potentialMember.addData("Ew0dL", Ew0dL, overwrite);

    potentialMember.addData("concAdd", concAdd, overwrite);
    potentialMember.addData("shadowing", shadowingInt, overwrite);
    potentialMember.addData("ConcEquilReset", resetInt, overwrite);
    potentialMember.addData("useConcEquil", concEquilInt, overwrite);
    potentialMember.addData("usediffusion", useDiffusionInt, overwrite);
    potentialMember.addData("useisotropicdiffusion", isotropicDiffusionInt, overwrite);
    potentialMember.addData("size", size.value(), overwrite);
    potentialMember.addData("heightmap", heightmap, overwrite);
    potentialMember.addData("ignisData", solver->mainLattice()->storedEventValues(), overwrite);
    potentialMember.addData("ignisEventDescriptions", solver->mainLattice()->outputEventDescriptions(), overwrite);
    if (acf)
    {
        potentialMember.addData("AutoCorr", autocorr.acf(), overwrite);
    }

    if (concEquil && !reset)
    {
        potentialMember.addData("muEq", cc.averageMu(), overwrite);
        potentialMember.addData("muEqError", cc.error(), overwrite);
    }

    return 0;

}

