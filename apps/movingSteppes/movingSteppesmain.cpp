#include <kMC>
#include <libconfig_utils/libconfig_utils.h>

using namespace libconfig;
using namespace kMC;


void initializeMovingSteppes(KMCSolver * solver, const Setting & root);

int main()
{

    Config cfg;
    wall_clock t;


    cfg.readFile("infiles/movingSteppes.cfg");

    const Setting & root = cfg.getRoot();


    KMCDebugger_SetFilename("movingSteppes");

    KMCDebugger_SetEnabledTo(getSetting<int>(root, "buildTrace") == 0 ? false : true);


    KMCSolver* solver = new KMCSolver(root);

    initializeMovingSteppes(solver, root);


    t.tic();

    solver->mainloop();

    cout << "Simulation ended after " << t.toc() << " seconds" << endl;


    KMCDebugger_DumpFullTrace();

    return 0;

}

class MovingSteppesEvent : public KMCEvent
{
public:

    MovingSteppesEvent() : KMCEvent("movingSteppes", "", true, true) {}

protected:

    void execute()
    {
        setValue(solver()->solverEvent()->totalTime());
    }

};

void initializeMovingSteppes(KMCSolver *solver, const Setting &root)
{

    solver->initializeFromXYZ("/home/jorgen/code/build-kMC-Desktop_Qt_5_2_1_GCC_64bit-Release/apps/movingSteppes/outfiles", 1659);
    return;

    const Setting & initCFG = getSetting(root, "Initialization");

    const uint & NX = solver->NX();
    const uint & NY = solver->NY();

    const uint & height  = getSetting<uint>(initCFG, "height");
    const uint & padding = getSetting<uint>(initCFG, "padding");

    for (uint z = 0; z < padding; ++z)
    {
        for (uint y = 0; y < NY; ++y)
        {
            for (uint x = 0; x < NX; ++x)
            {
                solver->forceSpawnParticle(x, y, z);
            }
        }
    }

    const uint layerwidth = NX/(2*height+1);

    for (uint z = 0; z < height; ++z)
    {
        for (uint y = 0; y < NY; ++y)
        {
            for (uint x = (NX + layerwidth)/2 + layerwidth*z; x < NX; ++x)
            {
                solver->forceSpawnParticle(x, y, padding + z);
            }

            for (uint x = 0; x < (NX-layerwidth)/2 - z*layerwidth; ++x)
            {
                solver->forceSpawnParticle(x, y, padding + z);
            }
        }
    }

//    solver->addEvent(new MovingSteppesEvent());

}
