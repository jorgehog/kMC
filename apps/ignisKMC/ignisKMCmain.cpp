#include <kMC>
#include <libconfig_utils/libconfig_utils.h>
#include <DCViz.h>

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


//    DCViz viz("/tmp/ignisEventsOut.arma");
//    viz.launch(true, 0.2, 30, 16);

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

    tempChange(const double T1, uint therm = 1) :
        KMCEvent("tempChange"),
        therm(therm),
        T1(T1)
    {

    }

    void initialize()
    {
        T0 = DiffusionReaction::beta();
        dT = (T1 - T0)/((eventLength/(double)therm - 1));

    }

    void execute()
    {
        if (nTimesExecuted%therm == 0)
        {
            DiffusionReaction::setBeta(T0 + dT*(nTimesExecuted/therm));
        }
    }



private:

    uint therm;

    double T0;
    const double T1;

    double dT;


};

class MeasureTemp : public KMCEvent
{
public:

    MeasureTemp() : KMCEvent("Temperature", "T*", true, true) {}

protected:

    void execute()
    {
        setValue(DiffusionReaction::beta());
    }

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

        solver()->forEachParticleDo([&cN, &c] (SoluteParticle *particle)
        {
            cN += particle->nNeighbors();
            c++;
        });

        setValue(cN/(double(c)));
    }

};

class TotalEnergy : public KMCEvent
{
public:

    TotalEnergy() : KMCEvent("TotalEnergy", "E*", true, true) {}

protected:

    void execute()
    {
        setValue(SoluteParticle::totalEnergy());
    }

};


class Debug : public KMCEvent
{
public:

    Debug() : KMCEvent("debug", "", true) {}

protected:

    void execute()
    {
        setValue(particles(0, 0));
    }
};



void initialize_ignisKMC(KMCSolver * solver, const Setting & root)
{

    const Setting & initCFG = getSurfaceSetting(root, "Initialization");

    const uint & NX = solver->NX();
    const uint & NY = solver->NY();
    const uint & NZ = solver->NZ();

    (void) NX;
    (void) NY;
    (void) NZ;

    const uint therm     = getSurfaceSetting<uint>(initCFG, "therm");
    const double endBeta = getSurfaceSetting<double>(initCFG, "endBeta");


    KMCEvent *heating = new tempChange(endBeta, therm);
    heating->setOffsetTime(MainLattice::nCycles/2 - 1);

    KMCEvent *cooling = new tempChange(DiffusionReaction::beta(), therm);
    cooling->setOnsetTime(MainLattice::nCycles/2);


    solver->addEvent(*heating);
    solver->addEvent(*cooling);

    solver->addEvent(new MeasureTemp());
    solver->addEvent(new AverageNeighbors());
    solver->addEvent(new TotalEnergy());


    solver->initializeSolutionBath();


}
