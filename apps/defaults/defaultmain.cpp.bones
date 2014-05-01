#include <kMC>
#include <libconfig_utils/libconfig_utils.h>

using namespace libconfig;
using namespace kMC;


void initialize__cname__(KMCSolver * solver, const Setting & root);

int main()
{

    Config cfg;
    wall_clock t;


    cfg.readFile("infiles/__name__.cfg");

    const Setting & root = cfg.getRoot();


    KMCDebugger_SetFilename("__name__");

    KMCDebugger_SetEnabledTo(getSetting<int>(root, "buildTrace") == 0 ? false : true);


    KMCSolver* solver = new KMCSolver(root);

    initialize__cname__(solver, root);


    t.tic();

    solver->mainloop();

    cout << "Simulation ended after " << t.toc() << " seconds" << endl;


    KMCDebugger_DumpFullTrace();

    return 0;

}

class __cname__Event : public KMCEvent
{
public:

    __cname__Event() : KMCEvent("__name__", "", true, true) {}

protected:

    void execute()
    {
        setValue(0);
    }

};

void initialize__cname__(KMCSolver *solver, const Setting &root)
{

    const Setting & initCFG = getSetting(root, "Initialization");

    const uint & NX = solver->NX();
    const uint & NY = solver->NY();
    const uint & NZ = solver->NZ();

    const uint & someValue = getSetting<uint>(initCFG, "someValue");

    (void)NX;
    (void)NY;
    (void)NZ;

    (void)someValue;

    solver->addEvent(new __cname__Event());

}
