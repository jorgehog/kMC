#include <kMC>
#include <libconfig_utils/libconfig_utils.h>

#include "WLMCEvent.h"

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


void initializeWLMC(KMCSolver *solver, const Setting &root)
{

    const Setting &initCFG = getSetting(root, "Initialization");

    const uint &nbins = getSetting<uint>(initCFG, "nbins");

    const uint &movesPerSampling = getSetting<uint>(initCFG, "movesPerSampling");

    const double &flatnessCriterion = getSetting<double>(initCFG, "flatnessCriterion");

    const uint &nbinsOverOverlap = getSetting<uint>(initCFG, "nbinsOverOverlap");
    const uint overlap = nbins/nbinsOverOverlap;

    const uint &nbinsOverMinWindowSize = getSetting<uint>(initCFG, "nbinsOverMinWindowSize");
    const uint minWindowSize = nbins/nbinsOverMinWindowSize;

    const uint &windowIncrementSize = getSetting<uint>(initCFG, "windowIncrementSize");

    const double &fStart = getSetting<double>(initCFG, "fStart");

    const double &fFinalMinusOne = getSetting<double>(initCFG, "fFinalMinusOne");
    const double fFinal = fFinalMinusOne + 1;


    KMCEvent *wlmc = new WLMCEvent(nbins,
                                   movesPerSampling,
                                   flatnessCriterion,
                                   overlap,
                                   minWindowSize,
                                   windowIncrementSize,
                                   fStart,
                                   fFinal);

    solver->swapMainSolverEventWith(wlmc);

}
