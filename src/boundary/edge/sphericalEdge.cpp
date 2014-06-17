#include "sphericalEdge.h"

#include "../../kmcsolver.h"
#include "../../site.h"

using namespace kMC;

SphericalEdge::SphericalEdge(const uint dimension, const uint orientation) :
    Boundary(dimension, orientation, Boundary::SphericalEdge)
{

}

void SphericalEdge::update()
{

}

void SphericalEdge::initialize()
{
    m_interface.reset();

    applyBoundaryTransform(NX(), NY(), NZ(), [this] (const uint &span, const uint &L, const uint &W)
    {
        m_interface.set_size(L, W);

        int W2 = W/2;
        int L2 = L/2;
        int S2 = span/2;

        double fl = 1.0/((L2-1)*(L2-1));
        double fw = 1.0/((W2-1)*(W2-1));

        for (int dl = -L2; dl < L2; ++dl)
        {
            for (int dw = -W2; dw < W2; ++dw)
            {
                uint ds = round((S2-1)*sqrt(1 - fl*dl*dl - fw*dw*dw));

                uint s = S2 + orientationAsSign()*ds;
                uint l = dl + L2;
                uint w = dw + W2;

                m_interface(l, w) = s;

                applyInverseBoundaryTransform(s, l, w, [this] (uint x, uint y, uint z)
                {
                    if (solver()->getSite(x, y, z)->isActive())
                    {
                        return;
                    }

                    solver()->forceSpawnParticle(x, y, z);
                });

            }
        }
    });

}
