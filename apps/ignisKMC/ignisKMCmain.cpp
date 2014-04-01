#include <kMC_ignis>
#include <libconfig_utils/libconfig_utils.h>

using namespace libconfig;
using namespace kMC;
using namespace ignis;


MainMesh *initialize_ignisKMC(KMCSolver * solver, const Setting & root);

int main()
{

    Config cfg;
    wall_clock t;


    cfg.readFile("infiles/ignisKMC.cfg");

    const Setting & root = cfg.getRoot();


    KMCDebugger_SetFilename("ignisKMC");

    KMCDebugger_SetEnabledTo(getSurfaceSetting<int>(root, "buildTrace") == 0 ? false : true);


    KMCSolver* solver = new KMCSolver(root);

    MainMesh* mainMesh = initialize_ignisKMC(solver, root);


    t.tic();

    mainMesh->eventLoop(solver->nCycles());

    cout << "Simulation ended after " << t.toc() << " seconds" << endl;


    KMCDebugger_DumpFullTrace();

    delete solver;


    return 0;

}


MainMesh *initialize_ignisKMC(KMCSolver * solver, const Setting & root)
{

    const Setting & initCFG = getSurfaceSetting(root, "Initialization");

    const uint & NX = solver->NX();
    const uint & NY = solver->NY();
    const uint & NZ = solver->NZ();

    mat topology;
    topology << 0 << NX << endr << 0 << NY << endr << 0 << NZ;

    solver->initializeSolutionBath();

    Particles *ens = new Particles(vec({0}));

    MainMesh *mainMesh = new MainMesh(topology, *ens);


    ReportProgress* a = new ReportProgress;
    mainMesh->addEvent(*a);

    SolverEvent *se = new SolverEvent(solver);

    mainMesh->addEvent(*se);


    return mainMesh;


}
