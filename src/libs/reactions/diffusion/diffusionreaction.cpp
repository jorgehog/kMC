#include "diffusionreaction.h"
#include "../../kmcsolver.h"

DiffusionReaction::DiffusionReaction(Site* destination) :
    Reaction(),
    destination(destination)
{

}

void DiffusionReaction::calcRate()
{
    uint ns = 0;
    uint nns = 0;

    uint xs = ((x() + x1())%NX)/2;
    uint ys = ((y() + y1())%NY)/2;
    uint zs = ((z() + z1())%NZ)/2;

    for (uint is = 0; is < 6; ++is) {

        uint I = (x()-2 + is + NX)%NX;
        for (uint js = 0; js < 6; ++js) {

            uint J = (y() - 2 + js + NY)%NY;
            for (uint ks = 0; ks < 6; ++ks) {

                uint K = (z() - 2 + ks + NZ)%NZ;

                double dx = (I - xs);
                double dy = (J - ys);
                double dz = (K - zs);

                double l2 = dx*dx + dy*dy + dz*dz;

                if ((l2 >= 1.25) && (l2 < 1.5)) {
                    ns++;
                } else if ((l2 >= 1.5) && (l2 < 3)) {
                    nns++;
                }

            }
        }

    }

    double Esp = nns*EspN + nns*EspNN;

    m_rate = mu*exp(-(reactionSite->E-Esp)/temperature);
}

bool DiffusionReaction::isActive()
{

    //Diffusion is active if the destination is empty
    return !destination->active();
}

void DiffusionReaction::execute()
{
    mainSolver->deactivateSite(reactionSite);
    mainSolver->activateSite(destination);
}
