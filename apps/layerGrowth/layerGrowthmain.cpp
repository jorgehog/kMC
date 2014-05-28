#include <kMC>
#include <libconfig_utils/libconfig_utils.h>

using namespace libconfig;
using namespace kMC;


void initialize_layerGrowth(KMCSolver * solver, const Setting & root);

int main()
{

    Config cfg;
    wall_clock t;


    cfg.readFile("infiles/layerGrowth.cfg");

    const Setting & root = cfg.getRoot();


    KMCDebugger_SetFilename("layerGrowth");

    KMCDebugger_SetEnabledTo(getSetting<int>(root, "buildTrace") == 0 ? false : true);


    KMCSolver* solver = new KMCSolver(root);

    initialize_layerGrowth(solver, root);


    t.tic();

    solver->mainloop();

    cout << "Simulation ended after " << t.toc() << " seconds" << endl;


    KMCDebugger_DumpFullTrace();

    return 0;

}

class LayerSize : public KMCEvent
{
public:

    LayerSize() : KMCEvent("LayerSize", "", true, true) {}

    void initialize()
    {
        m_z = 0;
        getLayerCount();
    }

    const uint &z() const
    {
        return m_z;
    }

    const uint &count() const
    {
        return m_count;
    }

protected:

    void execute()
    {
        if (lastReaction() == NULL)
        {
            return;
        }

        while(isFullLayer() && m_z != solver()->NZ())
        {
            m_z++;
            getLayerCount();
        }

        if (prevZ() == m_z)
        {
            if (newZ() != m_z)
            {
                m_count--;
            }
        }

        else if (newZ() == m_z)
        {
            if (prevZ() != m_z)
            {
                m_count++;
            }
        }

        setValue(m_count/(solver()->NX()*(double)solver()->NY()));
    }

private:

    void getLayerCount()
    {
        m_count = 0;
        for (uint x = 0; x < solver()->NX(); ++x)
        {
            for (uint y = 0; y < solver()->NY(); ++y)
            {
                if (solver()->getSite(x, y, m_z)->isActive())
                {
                    m_count++;
                }
            }
        }
    }

    bool isFullLayer()
    {
        return m_count == solver()->NX()*solver()->NY();
    }

    uint newZ() const
    {
        return lastReaction()->reactant()->z();
    }

    uint prevZ() const
    {
        return newZ() - lastReaction()->path(2);
    }

    uint m_z;
    uint m_count;

};

class ClusterNess : public KMCEvent
{
public:

    ClusterNess(const LayerSize *layerEvent) :
        KMCEvent("Clusterness", "", true, true),
        m_layerEvent(layerEvent)
    {

    }

protected:

    void execute()
    {
        uint clusterness = 0;
        uint N = 0;

        Site *currentSite;
        SoluteParticle *particle;

        for (uint x = 0; x < NX(); ++x)
        {
            for (uint y = 0; y < NY(); ++y)
            {
                currentSite = solver()->getSite(x, y, m_layerEvent->z());

                if (currentSite->isActive())
                {
                    N++;
                    particle = currentSite->associatedParticle();

                    for (int dx = -1; dx <= 1; ++dx)
                    {
                        for (int dy = -1; dy <= 1; ++dy)
                        {
                            if (dx == dy && dy == 0)
                            {
                                continue;
                            }

                            if (particle->neighborhood(dx, dy, 0)->isActive())
                            {
                                clusterness++;
                            }
                        }
                    }
                }
            }
        }

        setValue(clusterness/(double)N/8);

    }

private:

    const LayerSize *m_layerEvent;

};

class RandomInsertion : public KMCEvent
{
public:

    RandomInsertion() : KMCEvent() {}

    void initialize()
    {
        m_nPrev = SoluteParticle::nSolutionParticles();
    }

protected:

    void execute()
    {

        if (m_nTimesExecuted%10 == 0)
        {
            if (SoluteParticle::nSolutionParticles() == 0)
            {
                solver()->insertRandomParticle();
            }
        }

        m_nPrev = SoluteParticle::nSolutionParticles();
    }

private:

    uint m_nPrev;

};

void initialize_layerGrowth(KMCSolver * solver, const Setting & root)
{

    const Setting & initCFG = getSetting(root, "Initialization");

    const uint & NX = solver->NX();
    const uint & NY = solver->NY();

    const uint nInitialLayers = getSetting<uint>(initCFG, "nInitialLayers");
    const double lconc = getSetting<double>(initCFG, "surfaceConcentration");



    //fill first n layers
    for (uint z = 0; z < nInitialLayers; ++z)
    {
        for (uint x = 0; x < NX; ++x)
        {
            for (uint y = 0; y < NY; ++y)
            {
                solver->forceSpawnParticle(x, y, z);
            }
        }
    }

    const uint nSurfaceParticles = lconc*NX*NY;

    const double deltaX = sqrt(1.0/lconc*NX/(double)NY);
    const double deltaY = deltaX*NY/(double)NX;

    uint x = 0;
    uint y = 0;

    uint N = 0;

    uint i = 0;
    uint j = 0;

    while (N != nSurfaceParticles)
    {
        solver->forceSpawnParticle(x, y, nInitialLayers);

        N++;
        i++;

        x = i*deltaX;

        if (x >= NX)
        {
            j++;
            y = j*deltaY;

            i = 0;
            x = 0;

            if (y > NY)
            {
                break;
            }
        }
    }

//    LayerSize *layerSize = new LayerSize();
//    solver->addEvent(layerSize);
//    solver->addEvent(new ClusterNess(layerSize));

//    solver->addEvent(new RandomInsertion());

}
