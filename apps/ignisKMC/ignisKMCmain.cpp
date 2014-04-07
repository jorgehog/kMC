#include <kMC>
#include <libconfig_utils/libconfig_utils.h>

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

    KMCDebugger_SetEnabledTo(getSurfaceSetting<int>(root, "buildTrace") == 0 ? false : true);


    KMCSolver* solver = new KMCSolver(root);

    initialize_ignisKMC(solver, root);


    t.tic();

    solver->mainloop();

    cout << "Simulation ended after " << t.toc() << " seconds" << endl;


    KMCDebugger_DumpFullTrace();

    delete solver;


    return 0;

}

class tempChange : public KMCEvent
{
public:

    tempChange(const double T1) :
        KMCEvent("tempChange", "T0", true, true),
        T1(T1)
    {

    }

    void initialize()
    {
        T0 = DiffusionReaction::beta();

        dT = (T1 - T0)/eventLength;

    }

    void execute()
    {

        DiffusionReaction::setBeta(T0 + dT*(nTimesExecuted+1));
        setValue(DiffusionReaction::beta()/T0);
    }



private:

    double T0;
    const double T1;

    double dT;


};

class AverageNeighbors : public KMCEvent
{
public:

    AverageNeighbors() : KMCEvent("avgN", "", true, true) {}

protected:

    void execute()
    {
        uint cN = 0;
        uint c  = 0;
        solver()->forEachActiveSiteDo([&cN, &c] (Site * currentSite)
        {
            cN += currentSite->nNeighbors();
            c++;
        });

        setValue(cN/(double(c)));
    }

};



void initialize_ignisKMC(KMCSolver * solver, const Setting & root)
{

    const Setting & initCFG = getSurfaceSetting(root, "Initialization");

    const uint & NX = solver->NX();
    const uint & NY = solver->NY();
    const uint & NZ = solver->NZ();

    KMCEvent *tChange       = new tempChange(1.5);
    KMCEvent *avgNeighbors  = new AverageNeighbors();

    solver->addEvent(*tChange);
    solver->addEvent(*avgNeighbors);
    solver->initializeSolutionBath();






}
