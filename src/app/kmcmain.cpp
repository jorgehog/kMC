
#include <kmcsolver.h>

#include <libconfig_utils/libconfig_utils.h>

using namespace libconfig;


int main()
{

    Config cfg;

    cfg.readFile("infiles/config.cfg");

    const Setting & root = cfg.getRoot();
    const Setting & test = getSetting(root, {"FirstLayer", "NextLayer"});

    int someVariable = getSurfaceSetting<int>(test, "someVariable");

    cout << someVariable << endl;
    return 0;
//    KMCSolver* solver = new KMCSolver(nx, ny, nz);
//    solver->run();



    return 0;
}
