
#include <kmcsolver.h>

#include <libconfig_utils/libconfig_utils.h>

using namespace libconfig;


int main()
{

//    Config cfg;

//    cfg.readFile("infiles/config.cfg");

//    const Setting & root = cfg.getRoot();

//    uint nx = getSetting<uint>(root, {"System", "NX"});
//    uint ny = getSetting<uint>(root, {"System", "NY"});
//    uint nz = getSetting<uint>(root, {"System", "NZ"});

    uint nx = 100;
    uint ny = 100;
    uint nz = 100;

    KMCSolver* solver = new KMCSolver(nx, ny, nz);
    solver->run();



    return 0;
}
