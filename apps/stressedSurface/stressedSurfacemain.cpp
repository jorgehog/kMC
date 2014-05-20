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

class SpawnParticle : public KMCEvent
{
public:

    SpawnParticle(uint nFree) :
        KMCEvent("spawnParticle", "", true, true),
        m_nFree(nFree)
    {

    }

protected:

    void execute()
    {
        uint n = SoluteParticle::nSolutionParticles();

        if (n > m_nFree)
        {
            while (SoluteParticle::nSolutionParticles() > m_nFree)
            {

                for (SoluteParticle *particle : solver()->particles())
                {
                    if (particle->isSolvant())
                    {
                        solver()->despawnParticle(particle);
                        break;
                    }
                }

            }

            solver()->getRateVariables();
        }
        else if (n < m_nFree)
        {
            while (SoluteParticle::nSolutionParticles() < m_nFree)
            {
                solver()->insertRandomParticle();
            }

            solver()->getRateVariables();
        }


        setValue(SoluteParticle::nSolutionParticles());
    }

private:

    const uint m_nFree;

};

void initializeStressedSurface(KMCSolver *solver, const Setting &root)
{

    const Setting &initCFG = getSetting(root, "Initialization");

    const uint &NZ = solver->NZ();

    const uint &nFreeWalkers = getSetting<uint>(initCFG, "nFreeWalkers");

    const uint &nEdgeLayers = getSetting<uint>(initCFG, "nEdgeLayers");

    const double &initialHeightRatio = getSetting<double>(initCFG, "initialHeightRatio");
    const uint height = (NZ - nEdgeLayers)*initialHeightRatio;


    const double &Es = getSetting<double>(initCFG, "Es");
    const double &r0 = getSetting<double>(initCFG, "r0");

    //strong interaction layer bottom
    solver->initializeLayers(nEdgeLayers, 0, 1, true);

    //initial crystal rim
    solver->initializeLayers(height, nEdgeLayers);

    //quick hack, will cleanup when I get time.
    SoluteParticle::ss = new InertWall(Site::boundaries(2, 1), Es, r0, DiffusionReaction::rPower(1, 1), DiffusionReaction::strength(1, 1), 0.1);

//    solver->addEvent(new SpawnParticle(nFreeWalkers));

}
