#include "quasidiffusion.h"

#include <libconfig_utils/libconfig_utils.h>

using namespace libconfig;


void initializeQuasi2DLoaded(KMCSolver * solver, const Setting & root);

int main()
{

    Config cfg;
    wall_clock t;


    cfg.readFile("infiles/Quasi2DLoaded.cfg");

    const Setting & root = cfg.getRoot();


    KMCDebugger_SetFilename("Quasi2DLoaded");

    KMCDebugger_SetEnabledTo(getSetting<int>(root, "buildTrace") == 0 ? false : true);

    KMCSolver::enableDumpLAMMPS(false);
    MainLattice::enableEventFile(false);

    KMCSolver* solver = new KMCSolver(root);

    initializeQuasi2DLoaded(solver, root);


    t.tic();

    solver->mainloop();

    cout << "Simulation ended after " << t.toc() << " seconds" << endl;


    KMCDebugger_DumpFullTrace();

    return 0;

}

class DumpHeighmap : public KMCEvent
{
public:

    DumpHeighmap(const ivec *heighmap) :
        KMCEvent(),
        m_heighmap(heighmap)
    {

    }

protected:

    void execute()
    {
        if (nTimesExecuted()%MainLattice::nCyclesPerOutput == 0)
        {
            m_heighmap->save(solver()->filePath() + "heighmap.arma");
        }
    }

private:

    const ivec* m_heighmap;

};


void initializeQuasi2DLoaded(KMCSolver *solver, const Setting &root)
{

    const Setting &initCFG = getSetting(root, "Initialization");

    const uint &NX = solver->NX();

    const uint &width = getSetting<uint>(initCFG, "width");

    ivec* heighmap = new ivec(NX);

    //Override standard diffusion. Necessary for quasi diffusive simulations.
    solver->setDiffusionType(KMCSolver::DiffusionTypes::None);


    for (uint site = 0; site < NX; ++site)
    {
        (*heighmap)(site) = KMC_RNG_UNIFORM()*width;

        SoluteParticle* particle = solver->forceSpawnParticle(site, 0, 0);
        particle->addReaction(new LeftHop(particle, heighmap));
        particle->addReaction(new RightHop(particle, heighmap));
        particle->addReaction(new Adsorbtion(particle, heighmap));
        particle->addReaction(new Desorbtion(particle, heighmap));

    }

    QuasiDiffusionReaction::initialize(heighmap->max());

    solver->addEvent(new DumpHeighmap(heighmap));

}
