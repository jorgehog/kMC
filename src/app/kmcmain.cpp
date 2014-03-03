
#include <kMC>

#include <libconfig_utils/libconfig_utils.h>

using namespace libconfig;


int main()
{

    Config cfg;

    cfg.readFile("infiles/knowncase.cfg");

    const Setting & root = cfg.getRoot();

    KMCDebugger_SetFilename("sameSaddleProblem");
    KMCSolver* solver = new KMCSolver(root);

    solver->run();



    return 0;
}
