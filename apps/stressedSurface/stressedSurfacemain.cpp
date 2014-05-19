#include <kMC>
#include <libconfig_utils/libconfig_utils.h>

using namespace libconfig;
using namespace kMC;


void initializeStressedSurface(KMCSolver * solver, const Setting & root);

int main()
{

    Config cfg;
    wall_clock t;


    cfg.readFile("infiles/stressedSurface.cfg");

    const Setting & root = cfg.getRoot();


    KMCDebugger_SetFilename("stressedSurface");

    KMCDebugger_SetEnabledTo(getSetting<int>(root, "buildTrace") == 0 ? false : true);


    KMCSolver* solver = new KMCSolver(root);

    initializeStressedSurface(solver, root);


    t.tic();

    solver->mainloop();

    cout << "Simulation ended after " << t.toc() << " seconds" << endl;


    KMCDebugger_DumpFullTrace();

    return 0;

}

class StressedSurfaceEvent : public KMCEvent
{
public:

    StressedSurfaceEvent() : KMCEvent("stressedSurface", "", true, true) {}

protected:

    void execute()
    {
        setValue(0);
    }

};

void initializeStressedSurface(KMCSolver *solver, const Setting &root)
{

    const Setting &initCFG = getSetting(root, "Initialization");

    const uint &NZ = solver->NZ();

    const double &initialHeightRatio = getSetting<double>(initCFG, "initialHeightRatio");
    const uint height = NZ*initialHeightRatio;

    const uint &nEdgeLayers = getSetting<uint>(initCFG, "nEdgeLayers");

    const double &Es = getSetting<double>(initCFG, "Es");
    const double &r0 = getSetting<double>(initCFG, "r0");


    solver->initializeLayers(height);

    solver->initializeLayers(nEdgeLayers, NZ - nEdgeLayers - 1, 1);


    StressedSurface *ss = new StressedSurface(Site::boundaries(2, 1), Es, r0, nEdgeLayers);
    SoluteParticle::ss = ss; //quick hack, will cleanup when I get time.


    solver->addEvent(new StressedSurfaceEvent());

}
