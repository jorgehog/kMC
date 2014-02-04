#include "site.h"
#include "reactions/reaction.h"
#include "reactions/diffusion/diffusionreaction.h"
#include "kmcsolver.h"


const uint &Site::nNeighborsLimit()
{
    return m_nNeighborsLimit;
}

const uint & Site::neighborhoodLength()
{
    return m_neighborhoodLength;
}


const ucube &Site::levelMatrix()
{
    return m_levelMatrix;
}


const ivec &Site::originTransformVector()
{
    return m_originTransformVector;
}


Site::Site(uint _x, uint _y, uint _z) :
    E(0),
    m_x(_x),
    m_y(_y),
    m_z(_z)
{

    m_nNeighbors.set_size(m_nNeighborsLimit);
    m_nNeighbors.zeros();

    neighborHood = new Site***[m_neighborhoodLength];

    for (uint i = 0; i < m_neighborhoodLength; ++i) {

        neighborHood[i] = new Site**[m_neighborhoodLength];

        for (uint j = 0; j < m_neighborhoodLength; ++j) {

            neighborHood[i][j] = new Site*[m_neighborhoodLength];

        }
    }


}

Site::~Site()
{
    m_activeReactions.clear();
    m_siteReactions.clear();
}

void Site::loadNeighborLimit(const Setting &setting)
{
    const uint  &limit = getSurfaceSetting<uint>(setting, "nNeighborsLimit");

    m_nNeighborsLimit = limit;
    m_neighborhoodLength = 2*m_nNeighborsLimit + 1;

    m_levelMatrix = zeros<ucube>(m_neighborhoodLength, m_neighborhoodLength, m_neighborhoodLength);

    m_originTransformVector = linspace<ivec>(-(int)m_nNeighborsLimit, m_nNeighborsLimit, m_neighborhoodLength);

    for (uint i = 0; i < m_neighborhoodLength; ++i)
    {

        for (uint j = 0; j < m_neighborhoodLength; ++j)
        {

            for (uint k = 0; k < m_neighborhoodLength; ++k)
            {

                if (i == m_nNeighborsLimit && j == m_nNeighborsLimit && k == m_nNeighborsLimit)
                {
                    continue;
                }

                m_levelMatrix(i, j, k) = getLevel(std::abs(m_originTransformVector(i)),
                                                  std::abs(m_originTransformVector(j)),
                                                  std::abs(m_originTransformVector(k)));
            }
        }
    }

}

void Site::addReaction(Reaction *reaction)
{
    reaction->setSite(this);

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

    for (Reaction* reaction : m_activeReactions) {
        reaction->calcRate();
    }
}

void Site::setSolverPtr(KMCSolver *solver)
{
    NX = solver->getNX();
    NY = solver->getNY();
    NZ = solver->getNZ();

    mainSolver = solver;
}

void Site::distanceTo(const Site *other, int &dx, int &dy, int &dz, bool absolutes) const
{


    dx = (other->x() + NX - m_x)%NX;
    dy = (other->y() + NY - m_y)%NY;
    dz = (other->z() + NZ - m_z)%NZ;

    if ((uint)abs(dx) > NX/2) {
        dx = -(int)(NX - dx);
    }

    if ((uint)abs(dy) > NY/2) {
        dy = -(int)(NY - dy);
    }

    if ((uint)abs(dz) > NZ/2) {
        dz = -(int)(NZ - dz);
    }

    if (absolutes) {
        dx = abs(dx);
        dy = abs(dy);
        dz = abs(dz);
    }

}

void Site::updateEnergy(Site *changedSite, int change)
{
    int xScaled;
    int yScaled;
    int zScaled;

    distanceTo(changedSite, xScaled, yScaled, zScaled, true);

    xScaled += Site::nNeighborsLimit();
    yScaled += Site::nNeighborsLimit();
    zScaled += Site::nNeighborsLimit();

    double dE = change*DiffusionReaction::potential()(xScaled, yScaled, zScaled);

    E += dE;

    m_totalEnergy += dE;

}


void Site::introduceNeighborhood()
{
    uint xTrans, yTrans, zTrans;
    for (uint i = 0; i < m_neighborhoodLength; ++i) {

        xTrans = (m_x + m_originTransformVector(i) + NX)%NX;

        for (uint j = 0; j < m_neighborhoodLength; ++j) {

            yTrans = (m_y + m_originTransformVector(j) + NY)%NY;

            for (uint k = 0; k < m_neighborhoodLength; ++k) {

                zTrans = (m_z + m_originTransformVector(k) + NZ)%NZ;

                neighborHood[i][j][k] = mainSolver->getSites()[xTrans][yTrans][zTrans];

            }
        }
    }
}

void Site::informNeighborhoodOnChange(int change)
{

    Site *neighbor;
    uint level;

    for (uint i = 0; i < m_neighborhoodLength; ++i) {

        for (uint j = 0; j < m_neighborhoodLength; ++j) {

            for (uint k = 0; k < m_neighborhoodLength; ++k) {


                neighbor = neighborHood[i][j][k];

                if (neighbor == this) {
                    assert(i == j && j == k && k == m_nNeighborsLimit);
                    continue;
                }

                level = m_levelMatrix(i, j, k);
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

    for (uint i = 0; i < m_neighborhoodLength; ++i) {

        for (uint j = 0; j < m_neighborhoodLength; ++j) {

            for (uint k = 0; k < m_neighborhoodLength; ++k) {

                neighbor = neighborHood[i][j][k];

                if (neighbor == this) {
                    assert(i == m_nNeighborsLimit && j == m_nNeighborsLimit && k == m_nNeighborsLimit);
                    continue;
                }


                if (neighbor->active()) {

                    level = m_levelMatrix(i, j, k);
                    m_nNeighbors(level)++;

                    E += DiffusionReaction::potential()(i, j, k);
                    m_totalEnergy += DiffusionReaction::potential()(i, j, k);

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


KMCSolver *Site::mainSolver;

uint Site::NX;
uint Site::NY;
uint Site::NZ;

uint Site::m_nNeighborsLimit;

uint Site::m_neighborhoodLength;

ucube Site::m_levelMatrix;
ivec Site::m_originTransformVector;

uint Site::m_totalActiveSites = 0;

double Site::m_totalEnergy = 0;
