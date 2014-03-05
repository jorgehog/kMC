
#include <kMC>

#include <libconfig_utils/libconfig_utils.h>

using namespace libconfig;


int main()
{

    Config cfg;

    cfg.readFile("infiles/config.cfg");

    const Setting & root = cfg.getRoot();

    KMCDebugger_SetFilename("sameSaddleProblem");
    KMCSolver* solver = new KMCSolver(root);

    wall_clock t;
    t.tic();
    solver->run();
    cout << "Simulation ended after " << t.toc() << " seconds" << endl;


    return 0;
}
