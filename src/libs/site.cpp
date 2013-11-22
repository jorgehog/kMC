#include "site.h"
#include "reactions/reaction.h"
#include "kmcsolver.h"

Site::Site(uint _x, uint _y, uint _z, KMCSolver *solver) :
    mainSolver(solver),
    m_x(_x),
    m_y(_y),
    m_z(_z)
{

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

    E = En*nNeighbours() + Enn*nNextNeighbours();

    for (Reaction* reaction : m_activeReactions) {
        reaction->calcRate();
    }
}

uint Site::nNeighbours()
{
    return mainSolver->nNeighbours(m_x, m_y, m_z);
}

uint Site::nNextNeighbours()
{
    return  mainSolver->nNextNeighbours(m_x, m_y, m_z);
}
