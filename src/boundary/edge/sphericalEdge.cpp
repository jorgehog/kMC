#include "sphericalEdge.h"

#include "../../kmcsolver.h"
#include "../../site.h"
#include "../../soluteparticle.h"


using namespace kMC;

SphericalEdge::SphericalEdge(const uint dimension, const uint orientation) :
    Boundary(dimension, orientation, Boundary::SphericalEdge)
{

}

void SphericalEdge::update()
{
    removeOutsiders();

    performConcentrationBoundaryConditionStep();
}

void SphericalEdge::initialize()
{
    m_interface.reset();

    applyBoundaryTransform(NX(), NY(), NZ(), [this] (const uint &H, const uint &L, const uint &W)
    {
        m_interface.set_size(L, W);

        int W2 = W/2;
        int L2 = L/2;
        int H2 = H/2;

        double fl = 1.0/((L2-1)*(L2-1));
        double fw = 1.0/((W2-1)*(W2-1));

        for (int dl = -L2; dl < L2; ++dl)
        {
            for (int dw = -W2; dw < W2; ++dw)
            {
                uint dh = round((H2-1)*sqrt(1 - fl*dl*dl - fw*dw*dw));

                uint h = H2 + orientationAsSign()*dh;
                uint l = dl + L2;
                uint w = dw + W2;

                m_interface(l, w) = h;

            }
        }
    });

}

void SphericalEdge::removeOutsiders()
{
    bool outside;

    vector<SoluteParticle*> outsiders;

    for (SoluteParticle *particle : SoluteParticle::affectedParticles())
    {
        uint h, l, w;

        applyBoundaryTransform(particle->x(), particle->y(), particle->z(), [&h, &l, &w] (uint _h, uint _l, uint _w)
        {
            h = _h;
            l = _l;
            w = _w;
        });


        uint hInterface = m_interface(l, w);

        outside = orientation() == Near ? h < hInterface : h > hInterface;

        if (outside)
        {
            //solver::despawnparticle removes the particle from the list we are currently
            //traversing, so we need to store it for later removal.
            outsiders.push_back(particle);
        }
    }

    for (SoluteParticle *particle : outsiders)
    {
        solver()->despawnParticle(particle);
    }
}
