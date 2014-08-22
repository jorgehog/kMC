#include <kMC>
#include <commonkmcevents.h>

#include <libconfig_utils/libconfig_utils.h>
#include <DCViz/include/DCViz.h>

#include <armadillo>

using namespace libconfig;
using namespace kMC;
using namespace ignis;


void initialize_ignisKMC(KMCSolver * solver, const Setting & root);

int main()
{

    Config cfg;
    wall_clock t;


    cfg.readFile("infiles/ignisKMC.cfg");

    const Setting & root = cfg.getRoot();


    KMCDebugger_SetFilename("ignisKMC");

    KMCDebugger_SetEnabledTo(getSetting<int>(root, "buildTrace") == 0 ? false : true);


    KMCSolver* solver = new KMCSolver(root);

    initialize_ignisKMC(solver, root);

    DCViz viz("/tmp/ignisEventsOut.arma");
    viz.launch(true, 10, 30, 16);

    t.tic();

    solver->mainloop();

    cout << "Simulation ended after " << t.toc() << " seconds" << endl;


    KMCDebugger_DumpFullTrace();


    return 0;

}

void initialize_ignisKMC(KMCSolver * solver, const Setting & root)
{

    const Setting & initCFG = getSetting(root, "Initialization");

    const uint therm     = getSetting<uint>(initCFG, "therm");
    const double betaCoolMax = getSetting<double>(initCFG, "betaCoolMax");
    const double betaHeatMin = getSetting<double>(initCFG, "betaHeatMin");


    KMCEvent *cooling = new tempChange(betaCoolMax, therm);
    cooling->setOffsetTime(MainLattice::nCycles/2 - 1);

    KMCEvent *heating = new tempChange(betaHeatMin, therm);
    heating->setOnsetTime(MainLattice::nCycles/2);


    solver->addEvent(*cooling);
    solver->addEvent(*heating);

    solver->addEvent(new MeasureTemp());
    solver->addEvent(new AverageNeighbors());
    solver->addEvent(new TotalEnergy());
    solver->addEvent(new TotalTime());

    solver->initializeSolutionBath();


}
