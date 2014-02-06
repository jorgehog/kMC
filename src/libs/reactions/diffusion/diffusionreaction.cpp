#include "diffusionreaction.h"
#include "../../kmcsolver.h"

DiffusionReaction::DiffusionReaction(Site *destination) :
    Reaction(),
    destination(destination)
{

}

void DiffusionReaction::loadPotential(const Setting &setting)
{
    rPower = getSurfaceSetting<double>(setting, "rPower");
    scale  = getSurfaceSetting<double>(setting, "scale");

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

    saddleCutoff = 2*Site::nNeighborsLimit();

    m_saddleTransformVector = linspace<ivec>(-Site::nNeighborsLimit() + 1, Site::nNeighborsLimit(), saddleCutoff);

}

double DiffusionReaction::getSaddleEnergy()
{
    double xs = ((x() + x1())%NX)/2.0;
    double ys = ((y() + y1())%NY)/2.0;
    double zs = ((z() + z1())%NZ)/2.0;

    vector<Site*> neighborSet;

    std::set_intersection(reactionSite->allneighbors.begin(), reactionSite->allneighbors.end(),
                   destination->allneighbors.begin(), destination->allneighbors.end(),
                   std::back_inserter(neighborSet));


    double Esp = 0;

    int X, Y, Z;
    int X1, Y1, Z1;
    int d = (int)Site::nNeighborsLimit();

    for (Site* targetSite : neighborSet) {

        targetSite->distanceTo(reactionSite, X, Y, Z, true);
        targetSite->distanceTo(destination, X1, Y1, Z1, true);

        if ((X > d  && X1 > d)
                || (Y > d && Y1 > d)
                || (Z > d && Z1 > d))
        {
            cout << X << " " << Y << " " << Z << endl;
            cout << targetSite->x() << " is too far away from " << x() << "?" << endl;
            cout << targetSite->y() << " is too far away from " << y() << "?" << endl;
            cout << targetSite->z() << " is too far away from " << z() << "?" << endl;
            cout << "--------------------------------------" << endl;
            exit(1);
        }

        if (targetSite != reactionSite) {
            if (targetSite->active()) {
                double dx = fabs(xs - targetSite->x());
                double dy = fabs(ys - targetSite->y());
                double dz = fabs(zs - targetSite->z());

                if (dx > 3*Site::nNeighborsLimit()/2) {
                    dx = NX - dx;
                }

                if (dy > 3*Site::nNeighborsLimit()/2) {
                    dy = NY - dy;
                }

                if (dz > 3*Site::nNeighborsLimit()/2) {
                    dz = NZ - dz;
                }

                double r = sqrt(dx*dx + dy*dy + dz*dz);

                assert(r >= 1/2. && "Saddle point is atleast this distance from another site.");
                Esp += scale/pow(r, rPower);

            }
        }
    }

    lastUsedEsp = Esp;

    return Esp;

}

void DiffusionReaction::calcRate()
{
    lastUsedE = reactionSite->getEnergy();
    m_rate = mu*exp(-beta*(reactionSite->getEnergy()-getSaddleEnergy()));
}

bool DiffusionReaction::isActive()
{

    //if diffusion leads to increased potential energy we decline.
    return !destination->isBlocked() || destination->isSurface();
}

void DiffusionReaction::execute()
{
    reactionSite->deactivate();
    destination->activate();
}

double DiffusionReaction::rPower;
double DiffusionReaction::scale;

uint DiffusionReaction::saddleCutoff;

cube DiffusionReaction::m_potential;
ivec DiffusionReaction::m_saddleTransformVector;
