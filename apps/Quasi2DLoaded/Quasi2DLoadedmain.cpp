#include "quasidiffusion.h"
#include <commonkmcevents.h>

#include <libconfig_utils/libconfig_utils.h>
#include <DCViz/include/DCViz.h>

#include <HDF5Wrapper/include/hdf5wrapper.h>

using namespace libconfig;


ivec *initializeQuasi2DLoaded(KMCSolver * solver, const Setting & root, const uint l, const double beta);


class heightRMS : public KMCEvent
{
public:

    heightRMS(const ivec &heightmap) :
        KMCEvent("heightRMS", "l0", true, true),
        m_heightmap(heightmap),
        m_L(heightmap.size())
    {

    }


protected:

    void execute()
    {

        double meanHeight = sum(m_heightmap)/double(m_heightmap.size());

        double RMS = 0;

        for (uint i = 0; i < m_L; ++i)
        {
            RMS += (m_heightmap(i) - meanHeight)*(m_heightmap(i) - meanHeight);
        }

        RMS /= m_L;

        RMS = std::sqrt(RMS);

        setValue(RMS);

    }

private:

    const ivec &m_heightmap;
    const uint m_L;

};



int main()
{

    Config cfg;
    wall_clock t;


    cfg.readFile("infiles/Quasi2DLoaded.cfg");

    const Setting & root = cfg.getRoot();


    KMCDebugger_SetFilename("Quasi2DLoaded");

    KMCDebugger_SetEnabledTo(getSetting<int>(root, "buildTrace") == 0 ? false : true);

    KMCSolver::enableDumpLAMMPS(false);
    MainLattice::enableEventFile(true, "ignisQuasi2Dloaded.arma", 10000, 100);

    KMCSolver* solver = new KMCSolver(root);

    H5Wrapper::Root h5root("Quasi2D.h5");

    const Setting &initCFG = getSetting(root, "Initialization");

    const uint &lStart = getSetting<uint>(initCFG, "lStart");
    const uint &lEnd   = getSetting<uint>(initCFG, "lEnd");
    const uint &lStep  = getSetting<uint>(initCFG, "lStep");

    const double &tStart = getSetting<double>(initCFG, "tStart");
    const double &tEnd   = getSetting<double>(initCFG, "tEnd");
    const uint &nTemps  = getSetting<uint>(initCFG, "nTemps");

    vec temps = linspace(tStart, tEnd, nTemps);
    uint N = temps.size();

    t.tic();

    for (uint l = lStart; l <= lEnd; l += lStep)
    {
        for (uint i = 0; i < N; ++i)
        {
            double t = temps(i);

            cout << "Running l b = " << l << " " << t << endl;

            ivec* heightmap = initializeQuasi2DLoaded(solver, root, l, t);

            H5Wrapper::Member &sizeMember = h5root.addMember(l);

            stringstream s;
            s << QuasiDiffusionReaction::potentialString();

            H5Wrapper::Member &potentialMember = sizeMember.addMember(s.str());

            solver->mainloop();
            potentialMember.addData("heightmap", *heightmap);
            potentialMember.addData("ignisData", KMCEvent::eventMatrix());
            potentialMember.addData("ignisEventDescriptions", KMCEvent::outputEventDescriptions());

            potentialMember.file()->flush(H5F_SCOPE_GLOBAL);
            solver->reset();

            delete heightmap;
        }
    }

    cout << "Simulation ended after " << t.toc() << " seconds" << endl;

    KMCDebugger_DumpFullTrace();

    return 0;

}

class DumpHeighmap : public KMCEvent
{
public:

    DumpHeighmap(const ivec &heighmap) :
        KMCEvent(),
        m_heighmap(heighmap),
        m_filename(solver()->filePath() + "heighmap.arma"),
        m_viz(new DCViz(m_filename))
    {

    }

    ~DumpHeighmap()
    {
        delete m_viz;
    }

    void initialize()
    {
        //        m_viz->launch(true, 1, 20, 3);
    }

protected:

    void execute()
    {
        if (nTimesExecuted()%MainLattice::outputSpacing() == 0)
        {
            m_heighmap.save(m_filename);
        }
    }

private:

    const ivec &m_heighmap;

    const string m_filename;

    DCViz *m_viz;

};


ivec* initializeQuasi2DLoaded(KMCSolver *solver, const Setting &root, const uint l, const double beta)
{
    (void) root;

    solver->resetBoxSize(l, 1, 1);
    DiffusionReaction::setBeta(beta);

    ivec* heighmap = new ivec(l, fill::zeros);

    //Override standard diffusion. Necessary for quasi diffusive simulations.
    solver->setDiffusionType(KMCSolver::DiffusionTypes::None);


//    heighmap->at(0) = 0;
    for (uint site = 0; site < l; ++site)
    {
        //        if (site != l - 1)
        //        {
        //            (*heighmap)(site + 1) = (*heighmap)(site) + round(-1 + 2*KMC_RNG_UNIFORM());
        //        }

//        (*heighmap)(site) = -1 + 3*KMC_RNG_UNIFORM();

        SoluteParticle* particle = solver->forceSpawnParticle(site, 0, 0);
        particle->addReaction(new LeftHop(particle, *heighmap));
        particle->addReaction(new RightHop(particle, *heighmap));
        particle->addReaction(new Deposition(particle, *heighmap));
        //        particle->addReaction(new Desorbtion(particle, *heighmap));

    }

    QuasiDiffusionReaction::initialize();

    solver->addEvent(new DumpHeighmap(*heighmap));
    solver->addEvent(new TotalTime());
    solver->addEvent(new heightRMS(*heighmap));

    return heighmap;

}
