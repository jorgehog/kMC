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



void initialize_centerCrystal(KMCSolver * solver, const Setting & root)
{

    solver->initializeCrystal(getSetting<double>(root, {"Initialization", "RelativeSeedSize"}));

    solver->initializeSolutionBath();

}
