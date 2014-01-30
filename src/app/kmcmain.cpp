
#include <kmcsolver.h>

#include <libconfig_utils/libconfig_utils.h>

using namespace libconfig;


int main()
{

    Config cfg;

    cfg.readFile("infiles/config.cfg");

    const Setting & root = cfg.getRoot();

    uint nx = getSetting<uint>(root, {"System", "NX"});
    uint ny = getSetting<uint>(root, {"System", "NY"});
    uint nz = getSetting<uint>(root, {"System", "NZ"});

    uint nCycles = getSetting<uint>(root, {"Solver", "nCycles"});

    double beta = getSetting<double>(root, {"System", "beta"});
    double rPower = getSetting<double>(root, {"System", "rPower"});

    double saturationLevel = getSetting<double>(root, {"Initialization", "SaturationLevel"});
    int seedSize = getSetting<int>(root, {"Initialization", "SeedSize"});

    KMCSolver* solver = new KMCSolver(nCycles, nx, ny, nz);

    solver->run();



    return 0;
}
