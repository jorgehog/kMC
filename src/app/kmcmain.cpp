
#include <kmcsolver.h>

#include <libconfig_utils/libconfig_utils.h>

using namespace libconfig;


int main()
{

    Config cfg;

    cfg.readFile("infiles/config.cfg");

    const Setting & root = cfg.getRoot();

    KMCSolver* solver = new KMCSolver(root);

    solver->run();



    return 0;
}
