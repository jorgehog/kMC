#include <kMC>
#include <libconfig_utils/libconfig_utils.h>

using namespace libconfig;
using namespace kMC;


void initialize_surfaceGrowth(KMCSolver * solver, const Setting & root);

int main()
{

    Config cfg;
    wall_clock t;


    cfg.readFile("infiles/surfaceGrowth.cfg");

    const Setting & root = cfg.getRoot();


    KMCDebugger_SetFilename("surfaceGrowth");

    KMCDebugger_SetEnabledTo(getSetting<int>(root, "buildTrace") == 0 ? false : true);


    KMCSolver* solver = new KMCSolver(root);

    initialize_surfaceGrowth(solver, root);


    t.tic();

    solver->mainloop();

    cout << "Simulation ended after " << t.toc() << " seconds" << endl;


    KMCDebugger_DumpFullTrace();

    delete solver;


    return 0;

}

void initialize_surfaceGrowth(KMCSolver * solver, const Setting & root)
{

    const Setting & initCFG = getSetting(root, "Initialization");

    const uint & NX = solver->NX();
    const uint & NY = solver->NY();

    uint toothWidth   = getSetting<uint>(initCFG, "toothWidth");
    uint toothSpacing = getSetting<uint>(initCFG, "toothSpacing");

    uint toothHeight           =   toothWidth + 1;
    uint fullToothSize         = 2*toothWidth + 1;
    uint toothCenterSeparation = 2*toothWidth + toothSpacing;

    uint toothCenterX = 0;
    uint toothCenterY;

    uint X, Y;

    while (toothCenterX < NX)
    {

        toothCenterY = 0;

        while (toothCenterY < NY)
        {

            //Loop Z coordinate
            for (uint Z = 0; Z < toothHeight; ++Z)
            {

                //i, j = tooth coordiates relative to corners
                for (uint i = Z; i < fullToothSize - Z; ++i)
                {

                    //if relative coordinates are outside of the box, we continue.
                    if (toothCenterX + i < toothWidth)
                    {
                        continue;
                    }

                    //transforming to box coordinates
                    X = toothCenterX + i - toothWidth;

                    //if absolute coordiates are outside the box, we continue
                    if (X >= NX)
                    {
                        continue;
                    }

                    for(uint j = Z; j < fullToothSize - Z; ++j)
                    {

                        if (toothCenterY + j < toothWidth)
                        {
                            continue;
                        }

                        Y = toothCenterY + j - toothWidth;

                        if (Y >= NY)
                        {
                            continue;
                        }

                        solver->forceSpawnParticle(X, Y, Z);

                    }
                }

            }

            toothCenterY += toothCenterSeparation;

        }

        toothCenterX += toothCenterSeparation;

    }


    solver->initializeSolutionBath();

}
