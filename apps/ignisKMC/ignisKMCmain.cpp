#include <kMC>
#include <libconfig_utils/libconfig_utils.h>
#include <DCViz/include/DCViz.h>

#include <armadillo>

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

    KMCDebugger_SetEnabledTo(getSetting<int>(root, "buildTrace") == 0 ? false : true);


    KMCSolver* solver = new KMCSolver(root);

    initialize_ignisKMC(solver, root);


//    DCViz viz("/tmp/ignisEventsOut.arma");
//    viz.launch(true, 10, 30, 16);

    t.tic();

    solver->mainloop();

    cout << "Simulation ended after " << t.toc() << " seconds" << endl;


    KMCDebugger_DumpFullTrace();


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

        for (SoluteParticle *particle : solver()->particles())
        {
            cN += particle->nNeighbors();
            c++;
        }

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

class Sort : public KMCEvent
{
public:

    Sort() : KMCEvent() {}

protected:

    void execute()
    {
        if ((nTimesExecuted+1)%250000 == 0)
        {
            solver()->sortReactionsByRate();
        }
    }

};

class CPUTime : public KMCEvent
{
public:

    CPUTime() : KMCEvent("CPUTime", "s", true, true) {}

    void initialize()
    {
        clock.tic();
    }

protected:


    void execute()
    {
        setValue(clock.toc());
    }

    arma::wall_clock clock;


};


void initialize_ignisKMC(KMCSolver * solver, const Setting & root)
{

    const Setting & initCFG = getSetting(root, "Initialization");

    const uint & NX = solver->NX();
    const uint & NY = solver->NY();
    const uint & NZ = solver->NZ();

    (void) NX;
    (void) NY;
    (void) NZ;

    const uint therm     = getSetting<uint>(initCFG, "therm");
    const double betaCoolMax = getSetting<double>(initCFG, "betaCoolMax");
    const double betaHeatMin = getSetting<double>(initCFG, "betaHeatMin");


    KMCEvent *cooling = new tempChange(betaCoolMax, therm);
    cooling->setOffsetTime(MainLattice::nCycles/2 - 1);

    KMCEvent *heating = new tempChange(betaHeatMin, therm);
    heating->setOnsetTime(MainLattice::nCycles/2);


//    solver->addEvent(*cooling);
//    solver->addEvent(*heating);

    solver->addEvent(new MeasureTemp());
    solver->addEvent(new AverageNeighbors());
    solver->addEvent(new TotalEnergy());
//    solver->addEvent(new Sort());
//    solver->addEvent(new CPUTime());

    solver->initializeSolutionBath();

    lattice *sublattice = new lattice({0, 0, 0, NX/2, NY/2, NZ/2}, "sublattice");
    sublattice->addEvent(new countAtoms<uint>());

    solver->addSubLattice(sublattice);



}
