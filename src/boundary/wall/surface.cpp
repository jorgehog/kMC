#include "surface.h"

#include "../../kmcsolver.h"
#include "../../debugger/debugger.h"

using namespace kMC;

Surface::Surface(const uint dimension, const uint orientation) :
    Boundary(dimension, orientation, Boundary::Surface)
{

}

Surface::~Surface()
{

}



void Surface::initialize()
{

    setupBoundarySites();

    for (Site * boundarySite : boundarySites())
    {
        if (!boundarySite->isActive())
        {
            boundarySite->spawnAsFixedCrystal();
        }

        else
        {
            KMCDebugger_Assert(boundarySite->particleState(), ==, ParticleStates::fixedCrystal, "surface boundary attempted initialized with active particles on the surface.");
        }

    }

}

void Surface::finalize()
{
    for (Site * boundarySite : boundarySites())
    {
        if (boundarySite->isFixedCrystalSeed() && boundarySite->isActive())
        {
            boundarySite->deactivate();
        }
    }
}
