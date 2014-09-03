#include "quasidiffusion.h"
#include "quasidiffusionevents.h"

#include <commonkmcevents.h>


#include <libconfig_utils/libconfig_utils.h>

#include <HDF5Wrapper/include/hdf5wrapper.h>

using namespace libconfig;

class CustomSolverEvent : public KMCEvent
{
public:

    CustomSolverEvent() : KMCEvent("CustomSolverEvent")
    {

    }

    void initialize()
    {
        solver()->getRateVariables();

        m_totalTime = 0;
    }

    const double & totalTime() const
    {
        return m_totalTime;
    }

protected:

    void execute()
    {
        R = solver()->kTot()*KMC_RNG_UNIFORM();

        choice = solver()->getReactionChoice(R);

        m_selectedReaction = solver()->allPossibleReactions().at(choice);
        KMCDebugger_SetActiveReaction(m_selectedReaction);

        m_selectedReaction->execute();

        Site::updateBoundaries();

        solver()->getRateVariables();

        m_totalTime -= Reaction::linearRateScale()*std::log(KMC_RNG_UNIFORM())/solver()->kTot();

        //To counter buildup of roundoff errors
        if (m_nTimesExecuted % 10000 == 0)
        {
            solver()->remakeAccuAllRates();
        }

    }

private:

    double R;

    double m_totalTime;

    uint choice;

    Reaction * m_selectedReaction;

};


ivec *initializeQuasi2DLoaded(KMCSolver * solver, const Setting & root, const uint l, const double beta);

int main()
{

    Config cfg;
    wall_clock t;


    cfg.readFile("infiles/Quasi2DLoaded.cfg");

    const Setting & root = cfg.getRoot();


    KMCDebugger_SetFilename("Quasi2DLoaded");

    KMCDebugger_SetEnabledTo(getSetting<int>(root, "buildTrace") == 0 ? false : true);

    KMCSolver::enableLocalUpdating(false);
    KMCSolver::enableDumpLAMMPS(false);
    MainLattice::enableEventFile(true, "ignisQuasi2Dloaded.arma", 1, 1);

    KMCSolver* solver = new KMCSolver(root);

    H5Wrapper::Root h5root("Quasi2D.h5");

    const Setting &initCFG = getSetting(root, "Initialization");

    const uint &lStart = getSetting<uint>(initCFG, "lStart");
    const uint &lEnd   = getSetting<uint>(initCFG, "lEnd");
    const uint &lStep  = getSetting<uint>(initCFG, "lStep");

    const double &tStart = getSetting<double>(initCFG, "tStart");
    const double &tEnd   = getSetting<double>(initCFG, "tEnd");
    const uint &nTemps   = getSetting<uint>(initCFG, "nTemps");

    vec temps = linspace(tStart, tEnd, nTemps);
    uint N = temps.size();

    t.tic();

    for (uint l = lStart; l <= lEnd; l += lStep)
    {
        for (uint i = 0; i < N; ++i)
        {
            double t = temps(i);

            cout << "Running l b = " << l << " " << t << endl;

            ivec* heightmap = initializeQuasi2DLoaded(solver, initCFG, l, t);

            H5Wrapper::Member &sizeMember = h5root.addMember(l);

            H5Wrapper::Member &potentialMember = sizeMember.addMember(QuasiDiffusionReaction::potentialString());

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

ivec* initializeQuasi2DLoaded(KMCSolver *solver, const Setting &initCFG, const uint l, const double beta)
{

    const uint &h0 = getSetting<uint>(initCFG, "h0");
    const uint &nCells = getSetting<uint>(initCFG, "nCells");
    const double &concentrationFieldLength = getSetting<double>(initCFG, "concentrationFieldLength");

    const double &Eb = getSetting<double>(initCFG, "Eb");
    const double &EsMax = getSetting<double>(initCFG, "EsMax")*Eb;
    const double &EsInit = getSetting<double>(initCFG, "EsInit")*Eb;

    const double &chemicalPotentialDifference = getSetting<double>(initCFG, "chemicalPotentialDifference");
    const double &boundaryConcentration = getSetting<double>(initCFG, "boundaryConcentration");
    const double &diffusivity = getSetting<double>(initCFG, "diffusivity");

    solver->resetBoxSize(l, 1, 1);
    DiffusionReaction::setBeta(beta);

    ivec* heighmap = new ivec(l, fill::zeros);

    //Override standard diffusion. Necessary for quasi diffusive simulations.
    solver->setDiffusionType(KMCSolver::DiffusionTypes::None);

    //BAD PRATICE WITH POINTERS.. WILL FIX..
    ConcentrationControl *cc = new ConcentrationControl3D(boundaryConcentration, diffusivity, nCells, concentrationFieldLength);
    MovingWall *wallEvent = new MovingWall(h0, EsMax, EsInit, *heighmap, *cc);

    for (uint site = 0; site < l; ++site)
    {
        SoluteParticle* particle = solver->forceSpawnParticle(site, 0, 0);
        particle->addReaction(new LeftHopPressurized(particle, *heighmap, Eb, *wallEvent));
        particle->addReaction(new RightHopPressurized(particle, *heighmap, Eb, *wallEvent));
        particle->addReaction(new Deposition(particle, *heighmap, Eb, *wallEvent, chemicalPotentialDifference, *cc));
        particle->addReaction(new Dissolution(particle, *heighmap, Eb, *wallEvent, *cc));

    }

    solver->addEvent(wallEvent);
    solver->addEvent(new DumpHeighmap(*heighmap));
    solver->addEvent(new TotalTime());
    solver->addEvent(new heightRMS(*heighmap));

    return heighmap;

}
