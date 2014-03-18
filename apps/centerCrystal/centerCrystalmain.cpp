#include <kMC>
#include <libconfig_utils/libconfig_utils.h>

using namespace libconfig;
using namespace kMC;


void initialize_centerCrystal(KMCSolver * solver);

int main()
{

    Config cfg;
    wall_clock t;


    cfg.readFile("infiles/config.cfg");

    const Setting & root = cfg.getRoot();


    KMCDebugger_SetFilename("centerCrystal");

    KMCDebugger_SetEnabledTo(getSurfaceSetting<int>(root, "buildTrace") == 0 ? false : true);


    KMCSolver* solver = new KMCSolver(root);

    initialize_centerCrystal(solver);


    t.tic();

    solver->run();

    cout << "Simulation ended after " << t.toc() << " seconds" << endl;


    KMCDebugger_DumpFullTrace();

    delete solver;


    return 0;

}


void initialize_centerCrystal(KMCSolver * solver)
{

}
