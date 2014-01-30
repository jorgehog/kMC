#include "diffusionreaction.h"
#include "../../kmcsolver.h"

DiffusionReaction::DiffusionReaction(Site *destination) :
    Reaction(),
    destination(destination)
{

}

void DiffusionReaction::loadPotential(const Setting &setting)
{
    double rPower = getSurfaceSetting<double>(setting, "rPower");
    double scale  = getSurfaceSetting<double>(setting, "scale");

    m_potential.set_size(Site::neighborhoodLength(),
                         Site::neighborhoodLength(),
                         Site::neighborhoodLength());

    for (uint i = 0; i < Site::neighborhoodLength(); ++i)
    {

        for (uint j = 0; j < Site::neighborhoodLength(); ++j)
        {

            for (uint k = 0; k < Site::neighborhoodLength(); ++k)
            {

                if (i == Site::nNeighborsLimit() && j == Site::nNeighborsLimit() && k == Site::nNeighborsLimit())
                {
                    continue;
                }

                m_potential(i, j, k) = 1.0/pow(pow(Site::originTransformVector()(i), 2)
                                               + pow(Site::originTransformVector()(j), 2)
                                               + pow(Site::originTransformVector()(k), 2), rPower/2);
            }
        }
    }

    m_potential *= scale;
}

void DiffusionReaction::calcRate()
{
    uint ns = 0;
    uint nns = 0;

    uint xs = ((x() + x1())%NX)/2;
    uint ys = ((y() + y1())%NY)/2;
    uint zs = ((z() + z1())%NZ)/2;

    //? IS THIS RIGHT?
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

    m_rate = mu*exp(-beta*(reactionSite->getEnergy()-Esp));
}

bool DiffusionReaction::isActive()
{

    //Diffusion is active if the destination is empty
    return !destination->active();
}

void DiffusionReaction::execute()
{
    reactionSite->deactivate();
    destination->activate();
}

cube DiffusionReaction::m_potential;
