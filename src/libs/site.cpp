#include "site.h"
#include "reactions/reaction.h"
#include "reactions/diffusion/diffusionreaction.h"
#include "kmcsolver.h"

Site::Site(uint _x, uint _y, uint _z, KMCSolver *solver) :
    mainSolver(solver),
    E(0),
    m_x(_x),
    m_y(_y),
    m_z(_z)
{

    m_nNeighbors.set_size(nNeighborsLimit);

    neighborHood = new Site***[neighborhoodLength];

    for (uint i = 0; i < neighborhoodLength; ++i) {

        neighborHood[i] = new Site**[neighborhoodLength];

        for (uint j = 0; j < neighborhoodLength; ++j) {

            neighborHood[i][j] = new Site*[neighborhoodLength];

        }
    }


}

void Site::addReaction(Reaction *reaction)
{
    reaction->setSite(this);
    reaction->setMainsolver(mainSolver);

    m_siteReactions.push_back(reaction);
}

void Site::updateReactions()
{

    m_activeReactions.clear();

    if (!m_active) {
        return;
    }

    for (Reaction* reaction : m_siteReactions) {
        if (reaction->isActive()) {
            m_activeReactions.push_back(reaction);
        }
    }

}

void Site::calculateRates()
{

//    E = En*nNeighbors() + Enn*nNeighbors(1);

    for (Reaction* reaction : m_activeReactions) {
        reaction->calcRate();
    }
}

void Site::updateEnergy(Site *changedSite, int change)
{
    uint xScaled = (Site::nNeighborsLimit + abs((int)m_x - (int)changedSite->x()))%mainSolver->NX;
    uint yScaled = (Site::nNeighborsLimit + abs((int)m_y - (int)changedSite->y()))%mainSolver->NY;
    uint zScaled = (Site::nNeighborsLimit + abs((int)m_z - (int)changedSite->z()))%mainSolver->NZ;

    E += change*DiffusionReaction::weights(xScaled, yScaled, zScaled);

}


void Site::introduceNeighborhood()
{
    uint xTrans, yTrans, zTrans;
    for (uint i = 0; i < neighborhoodLength; ++i) {

        xTrans = (m_x + originTransformVector(i) + mainSolver->NX)%mainSolver->NX;

        for (uint j = 0; j < neighborhoodLength; ++j) {

            yTrans = (m_y + originTransformVector(j) + mainSolver->NY)%mainSolver->NY;

            for (uint k = 0; k < neighborhoodLength; ++k) {

                zTrans = (m_z + originTransformVector(k) + mainSolver->NZ)%mainSolver->NZ;

                neighborHood[i][j][k] = mainSolver->getSites()[xTrans][yTrans][zTrans];

            }
        }
    }
}

void Site::informNeighborhoodOnChange(int change)
{

    Site *neighbor;
    uint level;

    for (uint i = 0; i < neighborhoodLength; ++i) {

        for (uint j = 0; j < neighborhoodLength; ++j) {

            for (uint k = 0; k < neighborhoodLength; ++k) {


                neighbor = neighborHood[i][j][k];

                if (neighbor == this) {
                    continue;
                }

                level = levelMatrix(i, j, k);
                neighbor->m_nNeighbors(level)+=change;

                neighbor->updateEnergy(this, change);

                neighbor->updateReactions();
                neighbor->calculateRates();

            }
        }
    }

}

void Site::countNeighbors()
{

    m_nNeighbors.zeros();

    Site *neighbor;
    uint level;

    for (uint i = 0; i < neighborhoodLength; ++i) {

        for (uint j = 0; j < neighborhoodLength; ++j) {

            for (uint k = 0; k < neighborhoodLength; ++k) {

                neighbor = neighborHood[i][j][k];

                if (neighbor == this) {
                    continue;
                }


                if (neighbor->active()) {

                    level = levelMatrix(i, j, k);
                    m_nNeighbors(level)++;

                    E += DiffusionReaction::weights(i, j, k);

                }

            }
        }
    }
}

uint Site::getLevel(uint i, uint j, uint k)
{
    uint m = i;

    if (j > i) {
        m =  j;
    }

    if (k > m) {
        m = k;
    }

    return m - 1;
}

const uint Site::neighborhoodLength = 2*Site::nNeighborsLimit + 1;

ucube Site::levelMatrix = zeros<ucube>(Site::neighborhoodLength, Site::neighborhoodLength, Site::neighborhoodLength);
ivec Site::originTransformVector = linspace<ivec>(-(int)Site::nNeighborsLimit, Site::nNeighborsLimit, Site::neighborhoodLength);

uint Site::totalActiveSites = 0;
