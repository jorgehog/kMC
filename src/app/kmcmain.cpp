#include <kMC>

#include <libconfig_utils/libconfig_utils.h>

using namespace libconfig;
using namespace kMC;

int main()
{

    Config cfg;

    cfg.readFile("infiles/config.cfg");

    const Setting & root = cfg.getRoot();

    KMCSolver* solver = new KMCSolver(root);

    wall_clock t;

    KMCDebugger_SetFilename("kMCrun");
    KMCDebugger_SetEnabledTo(getSurfaceSetting<int>(root, "buildTrace") == 0 ? false : true);

    t.tic();
    solver->run();
    cout << "Simulation ended after " << t.toc() << " seconds" << endl;

    KMCDebugger_DumpFullTrace();
    delete solver;

    return 0;
}


//Introduce conditions for 'crystallizing' (and its inverse) as something easily replaced
//Same for diffusion. Introduce conditions for a blocked diffusion path (beyond weather site is occupied).
//Load into static  function as a lambda expression? Yes I think so!.
