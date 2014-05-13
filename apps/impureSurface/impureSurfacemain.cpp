#include <kMC>
#include <libconfig_utils/libconfig_utils.h>

using namespace libconfig;
using namespace kMC;


void initializeImpureSurface(KMCSolver * solver, const Setting & root);

int main()
{

    Config cfg;
    wall_clock t;


    cfg.readFile("infiles/impureSurface.cfg");

    const Setting & root = cfg.getRoot();


    KMCDebugger_SetFilename("impureSurface");

    KMCDebugger_SetEnabledTo(getSetting<int>(root, "buildTrace") == 0 ? false : true);

    KMCSolver::enableDumpLAMMPS(true);
    KMCSolver* solver = new KMCSolver(root);

    initializeImpureSurface(solver, root);


    t.tic();

    solver->mainloop();

    cout << "Simulation ended after " << t.toc() << " seconds" << endl;


    KMCDebugger_DumpFullTrace();

    return 0;

}

class ImpureSurfaceEvent : public KMCEvent
{
public:

    ImpureSurfaceEvent() : KMCEvent("impureSurface", "", true, true) {}

protected:

    void execute()
    {
        setValue(0);
    }

};

void initializeImpureSurface(KMCSolver *solver, const Setting &root)
{

    const Setting & initCFG = getSetting(root, "Initialization");

    const uint & NX = solver->NX();
    const uint & NY = solver->NY();

    const double & impurityDensity = getSetting<double>(initCFG, "impurityDensity");
    const uint & nInitialLayers = getSetting<uint>(initCFG, "nInitialLayers");

    const uint totalImpurities = impurityDensity*NX*NY;

    for (uint x = 0; x < NX; ++x)
    {
        for (uint y = 0; y < NY; ++y)
        {
            for (uint z = 0; z < nInitialLayers; ++z)
            {
                solver->forceSpawnParticle(x, y, z);
            }
        }
    }

    const double deltaX = sqrt(1.0/impurityDensity*NX/(double)NY);
    const double deltaY = deltaX*NY/(double)NX;

    uint x = 0;
    uint y = 0;

    uint N = 0;

    uint i = 0;
    uint j = 0;

    while (N != totalImpurities)
    {
        solver->forceSpawnParticle(x, y, nInitialLayers, 1);

        N++;
        i++;

        x = i*deltaX;

        if (x >= NX)
        {
            j++;
            y = j*deltaY;

            i = 0;
            x = 0;

            if (y >= NY)
            {
                break;
            }
        }
    }

    for (uint x = 0; x < NX; ++x)
    {
        for (uint y = 0; y < NY; ++y)
        {
            if (!solver->getSite(x, y, nInitialLayers)->isActive())
            {
                solver->forceSpawnParticle(x, y, nInitialLayers);
            }
        }
    }

    solver->addEvent(new ImpureSurfaceEvent());

}
