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
    uint xi = 0;

    if (orientation() == 1)
    {
        xi = span() - 1;
    }

    if (dimension() == X)
    {

        for (uint y = 0; y < NY(); ++y) {
            for (uint z = 0; z < NZ(); ++z) {
                if (!mainSolver()->getSite(xi, y, z)->isActive())
                {

                    mainSolver()->getSite(xi, y, z)->spawnAsFixedCrystal();
                }
            }
        }

    }

    else if (dimension() == Y)
    {

        for (uint x = 0; x < NX(); ++x) {
            for (uint z = 0; z < NZ(); ++z) {
                if (!mainSolver()->getSite(x, xi, z)->isActive())
                {
                    mainSolver()->getSite(x, xi, z)->spawnAsFixedCrystal();
                }
            }
        }

    }

    else
    {
        KMCDebugger_Assert(dimension(), ==, Z, "This else should always correspond to dim=2");

        for (uint x = 0; x < NX(); ++x) {
            for (uint y = 0; y < NY(); ++y) {
                if (!mainSolver()->getSite(x, y, xi)->isActive())
                {
                    mainSolver()->getSite(x, y, xi)->spawnAsFixedCrystal();
                }
            }
        }

    }


}
