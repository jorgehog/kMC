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

    SpawnParticle() :
        KMCEvent("spawnParticle", "", true, true)
    {

    }

    void initialize()
    {
        m_maxAbsFlux = 10;
    }

protected:

    void execute()
    {
        if (nTimesExecuted%100 != 0)
        {
            return;
        }

        uint N = std::round(solver()->targetConcentration()*SoluteParticle::getCurrentSolvantVolume()) + 0.00001;
        uint c = 0;

        if (SoluteParticle::nSolutionParticles() > N)
        {
            while (SoluteParticle::nSolutionParticles() > N && c < m_maxAbsFlux)
            {

                for (SoluteParticle *particle : solver()->particles())
                {
                    if (particle->isSolvant())
                    {
                        solver()->despawnParticle(particle);
                        c++;
                        break;
                    }
                }

            }

            solver()->getRateVariables();
        }
        else if (SoluteParticle::nSolutionParticles() < N)
        {
            while (SoluteParticle::nSolutionParticles() < N && c < m_maxAbsFlux)
            {
                solver()->insertRandomParticle();
                c++;
            }

            solver()->getRateVariables();
        }

        setValue(double(c)/m_maxAbsFlux);
    }

private:

    uint m_maxAbsFlux;

};

void initializeStressedSurface(KMCSolver *solver, const Setting &root)
{

    const Setting &initCFG = getSetting(root, "Initialization");

    const uint &NZ = solver->NZ();

    const double &initialHeightRatio = getSetting<double>(initCFG, "initialHeightRatio");
    const uint height = (NZ - 1)*initialHeightRatio;


    const double &Es = getSetting<double>(initCFG, "Es");
    const double &r0 = getSetting<double>(initCFG, "r0");

    //strong interaction layer bottom
    solver->initializeLayers(1, 0, 2, true);

    //initial crystal rim
    solver->initializeLayers(height, 1);

    //quick hack, will cleanup when I get time.
    SoluteParticle::ss = new InertWall(Site::boundaries(2, 1), Es, r0,
                                       DiffusionReaction::rPower(1, 1),
                                       DiffusionReaction::strength(1, 1), 0.1);

    solver->addEvent(new SpawnParticle());

    solver->initializeSolutionBath();

}
