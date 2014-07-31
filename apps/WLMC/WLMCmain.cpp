#include <kMC>
#include <BADAss/badass.h>
#include <libconfig_utils/libconfig_utils.h>

#include "kmcwlmcsystem.h"

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

    KMCDebugger_DumpFullTrace();

    return 0;

}


void initializeWLMC(KMCSolver *solver, const Setting &root)
{

    const Setting &initCFG = getSetting(root, "Initialization");

    const uint &adaptiveWindows = getSetting<uint>(initCFG, "adaptiveWindows");

    const uint &nParticlesStart = getSetting<uint>(initCFG, "nParticlesStart");
    const uint &nParticlesStop = getSetting<uint>(initCFG, "nParticlesStop");

    BADAss(nParticlesStart, <=, nParticlesStop);

    const uint &movesPerSampling = getSetting<uint>(initCFG, "movesPerSampling");

    const uint &nbins = getSetting<uint>(initCFG, "nbins");

    const uint &nbinsOverOverlap = getSetting<uint>(initCFG, "nbinsOverOverlap");
    const uint overlap = nbins/nbinsOverOverlap;

    const uint &nbinsOverMinWindowSize = getSetting<uint>(initCFG, "nbinsOverMinWindowSize");
    const uint minWindowSize = nbins/nbinsOverMinWindowSize;


    const double &flatnessCriterion = getSetting<double>(initCFG, "flatnessCriterion");

    const double &fStart = getSetting<double>(initCFG, "fStart");

    const double &fFinalMinusOne = getSetting<double>(initCFG, "fFinalMinusOne");
    const double fFinal = fFinalMinusOne + 1;


    KMCWLMCSystem *system;

    WLMC::Window *mainWindow;

    for (uint nParticles = nParticlesStart; nParticles <= nParticlesStop; ++nParticles)
    {
        solver->insertRandomParticles(nParticles);

        system = new KMCWLMCSystem(solver,
                                   movesPerSampling,
                                   flatnessCriterion,
                                   overlap,
                                   minWindowSize);

        mainWindow = system->execute(nbins, adaptiveWindows, fStart, fFinal);

        //do stuff.. setup folder structures etc. for loading later on.. use boost filesystem

        delete system;
        delete mainWindow;

        solver->clearParticles();
    }






}
