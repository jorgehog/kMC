#include "soluteparticle.h"

#include "kmcsolver.h"

#include "boundary/boundary.h"

#include "debugger/debugger.h"

#include "reactions/diffusion/diffusionreaction.h"


using namespace kMC;


SoluteParticle::SoluteParticle() :
    m_particleState(ParticleStates::solvant),
    m_x(UNSET_UINT),
    m_y(UNSET_UINT),
    m_z(UNSET_UINT),
    m_nNeighborsSum(0),
    m_energy(0),
    m_ID(ID_count++)
{

    initializeDiffusionReactions();

    setVectorSizes();


    refCounter++;

}

SoluteParticle::~SoluteParticle()
{

    KMCDebugger_Assert(refCounter, !=, 0);

    clearAllReactions();

    m_nNeighbors.clear();

    refCounter--;

}


void SoluteParticle::setNewPosition(const uint x, const uint y, const uint z)
{    

    KMCDebugger_Assert(solver()->particle(x, y, z), ==, NULL, "particle already present at site.", info());

    KMCDebugger_Assert(x, <, NX(), "mismatch in coordiantes. ", info());
    KMCDebugger_Assert(y, <, NY(), "mismatch in coordiantes. ", info());
    KMCDebugger_Assert(z, <, NZ(), "mismatch in coordiantes. ", info());

    m_x = x;
    m_y = y;
    m_z = z;

    m_solver->registerParticle(this);

    setupAllNeighbors();

    setNewParticleState(detectParticleState());

    forEachNeighborSiteDo_sendIndices([this] (SoluteParticle *neighbor, uint i, uint j, uint k)
    {
        neighbor->addNeighbor(this, i, j, k);
    });



    markAsAffected();

    forEachActiveReactionDo([] (Reaction *reaction)
    {
        reaction->forceUpdateFlag(Reaction::defaultUpdateFlag);
    });


    KMCDebugger_MarkPre("void");
    KMCDebugger_PushImplication(this, "enabled");
    KMCDebugger_MarkPartialStep("PARTICLE ACTIVATED");

}

void SoluteParticle::removeOldPosition()
{
    KMCDebugger_AssertBool(isActive(), "particle not present at site.", info());

    KMCDebugger_MarkPre(particleStateName());
    KMCDebugger_PushImplication(this, "disabled");

    m_solver->removeParticle(this);

    forEachNeighborSiteDo_sendIndices([this] (SoluteParticle *neighbor, uint i, uint j, uint k)
    {
        neighbor->removeNeighbor(this, i, j, k);
    });


    m_totalParticles(particleState())--;

    m_totalEnergy -= m_energy;

    m_energy = 0;


    KMCDebugger_MarkPartialStep("PARTICLE DISABLED");

}

bool SoluteParticle::isActive() const
{
    return solver()->isRegisteredParticle(this);
}

void SoluteParticle::setMainSolver(KMCSolver *solver)
{
    m_solver = solver;
}

void SoluteParticle::popAffectedParticle(SoluteParticle *particle)
{
    auto idx = m_affectedParticles.find(particle);

    if (idx == m_affectedParticles.end())
    {
        return;
    }

    KMCDebugger_AssertBool(particle->isAffected());
    KMCDebugger_PopAffected(particle);

    m_affectedParticles.erase(idx);
}

void SoluteParticle::updateAffectedParticles()
{

    KMCDebugger_PushTraces();

    for (SoluteParticle* particle : m_affectedParticles)
    {
        if (particle == NULL)
        {
            continue;
        }

        particle->updateReactions();
    }

    clearAffectedParticles();

}

double SoluteParticle::getCurrentSolvantVolume()
{
    return NX()*NY()*NZ() - nCrystals() - nSurfaces();
}





void SoluteParticle::setParticleState(int newState)
{

    m_totalParticles(m_particleState)--;
    m_totalParticles(newState)++;

    m_particleState = newState;

}


//All reactions must be legal if site is allowed to spawn.
bool SoluteParticle::isLegalToSpawn(const uint x, const uint y, const uint z) const
{
    for (int i = -1; i <= 1; ++i)
    {
        for (int j = -1; j <= 1; ++j)
        {
            for (int k = -1; k <= 1; ++k)
            {
                if (i == j && j == k && k == 0)
                {
                    continue;
                }

                else if (Site::neighborhood(x, y, z, i, j, k) == NULL)
                {
                    continue;
                }

                else if (Site::neighborhood(x, y, z, i, j, k)->isActive())
                {
                    return false;
                }
            }
        }

    }

    return true;
}


void SoluteParticle::clearAllReactions()
{

    for (Reaction * reaction : m_reactions)
    {
        delete reaction;
    }

    m_reactions.clear();

}



const string SoluteParticle::info(int xr, int yr, int zr, string desc) const
{
    stringstream s;

    s << "SoluteParticle" << m_ID << "@";
    s << Site::info(m_x, m_y, m_z, xr, yr, zr, desc);

    return s.str();

}

void SoluteParticle::setZeroEnergy()
{
    KMCDebugger_Assert(m_nNeighborsSum, ==, 0, "Energy is not zero.", info());
    m_energy = 0;

}

void SoluteParticle::setVectorSizes()
{

    m_nNeighbors.resize(Site::nNeighborsLimit());

}

uint SoluteParticle::nActiveReactions() const
{

    uint nActiveReactions = 0;

    forEachActiveReactionDo([&nActiveReactions] (Reaction * r)
    {
        (void) r;
        nActiveReactions++;
    });

    return nActiveReactions;
}

void SoluteParticle::markAsAffected()
{
    m_affectedParticles.insert(this);
}

void SoluteParticle::updateReactions()
{

    for (Reaction* reaction : m_reactions)
    {
        if (reaction->isAllowed())
        {
            reaction->calcRate();
            reaction->setLastUsedEnergy();
            reaction->resetUpdateFlag();
        }

        else
        {
            reaction->disable();
        }
    }
}

void SoluteParticle::setupAllNeighbors()
{

    KMCDebugger_Assert(m_energy, ==, 0, "Particle energy should be cleared prior to this call.", info());

    m_nNeighbors.zeros();

    m_nNeighborsSum = 0;

    forEachNeighborSiteDo_sendIndices([this] (SoluteParticle *neighbor, uint i, uint j, uint k)
    {
        this->addNeighbor(neighbor, i, j, k);
    });



}

void SoluteParticle::removeNeighbor(SoluteParticle *neighbor, uint i, uint j, uint k)
{
    _updateNeighborProps(-1, neighbor, i, j, k);
}

void SoluteParticle::addNeighbor(SoluteParticle *neighbor, uint i, uint j, uint k)
{
    _updateNeighborProps(+1, neighbor, i, j, k);
}

void SoluteParticle::_updateNeighborProps(const int sign, const SoluteParticle * neighbor, uint i, uint j, uint k)
{
    KMCDebugger_MarkPre(neighbor->particleStateName());

    const uint &level = Site::levelMatrix(i, j, k);

    m_nNeighbors(level) += sign;

    m_nNeighborsSum += sign;

    double dE = sign*DiffusionReaction::potential(i, j, k);

    m_energy += dE;

    m_totalEnergy += dE;

    if (level == 0)
    {
        changeParticleState(detectParticleState());
    }

    forEachActiveReactionDo([&level, &neighbor] (Reaction *reaction)
    {
        reaction->setDirectUpdateFlags(neighbor, level);
    });

    markAsAffected();

    KMCDebugger_PushImplication(neighbor, neighbor->particleStateName().c_str());
}

void SoluteParticle::distanceTo(const SoluteParticle *other, int &dx, int &dy, int &dz, bool absolutes) const
{
    Site::distanceBetween(m_x, m_y, m_z, other->x(), other->y(), other->z(), dx, dy, dz, absolutes);
}


void SoluteParticle::forEachActiveReactionDo(function<void (Reaction *)> applyFunction) const
{

    for (Reaction * reaction : m_reactions)
    {
        if (reaction->isAllowed())
        {
            applyFunction(reaction);
        }
    }
}

void SoluteParticle::forEachActiveReactionDo_sendIndex(function<void (Reaction *, uint)> applyFunction) const
{

    uint i = 0;

    for (Reaction * reaction : m_reactions)
    {
        if (reaction->isAllowed())
        {
            applyFunction(reaction, i);
        }

        i++;
    }
}


void SoluteParticle::initializeDiffusionReactions()
{


    KMCDebugger_Assert(m_reactions.size(), ==, 0, "Sitereactions are already set", info());

    //For each site, loop over all closest neighbors
    for (int dx = -1; dx <= 1; ++dx)
    {
        for (int dy = -1; dy <= 1; ++dy)
        {
            for (int dz = -1; dz <= 1; ++dz)
            {

                if (dx == 0 && dy == 0 && dz == 0)
                {
                    continue;
                }

                DiffusionReaction* diffusionReaction = new DiffusionReaction(this, dx, dy, dz);
                addReaction(diffusionReaction);

                m_diffusionReactions[dx + 1][dy + 1][dz + 1] = diffusionReaction;

            }
        }
    }
}

int SoluteParticle::detectParticleState()
{

    if (qualifiesAsCrystal())
    {
        return ParticleStates::crystal;
    }

    else if (isSolvant())
    {
        return ParticleStates::solvant;
    }

    KMCDebugger_AssertBool(qualifiesAsSurface());

    return ParticleStates::surface;

}


double SoluteParticle::potentialBetween(const SoluteParticle *other)
{
    int X, Y, Z;

    distanceTo(other, X, Y, Z, true);

    X += Site::nNeighborsLimit();
    Y += Site::nNeighborsLimit();
    Z += Site::nNeighborsLimit();

    return DiffusionReaction::potential(X, Y, Z);
}

uint SoluteParticle::maxDistanceTo(const SoluteParticle *other) const
{
    return Site::maxDistanceBetween(m_x, m_y, m_z, other->x(), other->y(), other->z());
}




void SoluteParticle::setZeroTotalEnergy()
{
    KMCDebugger_Assert(SoluteParticle::nParticles(), ==, 0);
    m_totalEnergy = 0;
}

uint SoluteParticle::nNeighborsSum() const
{
    KMCDebugger_Assert(sum(m_nNeighbors), ==, m_nNeighborsSum, "Should be identical.", info());
    return m_nNeighborsSum;
}

void SoluteParticle::clearAll()
{
    KMCDebugger_Assert(refCounter, ==, 0, "cannot clear static members with object alive.");

    m_totalParticles.zeros();

    setZeroTotalEnergy();

    clearAffectedParticles();

    ID_count = 0;
}

void SoluteParticle::clearAffectedParticles()
{
    m_affectedParticles.clear();
}


double SoluteParticle::getCurrentConcentration()
{
    if (nSolutionParticles() == 0)
    {
        return 0;
    }

    return static_cast<double>(nSolutionParticles())/getCurrentSolvantVolume();
}

double SoluteParticle::getCurrentRelativeCrystalOccupancy()
{
    return static_cast<double>(nCrystals() + nSurfaces())/(NX()*NY()*NZ());
}

void SoluteParticle::setNewParticleState(int newState)
{

    KMCDebugger_MarkPre(particleStateName());

    m_totalParticles(newState)++;


    m_particleState = newState;

    KMCDebugger_PushImplication(this, particleStateName().c_str());

}


void SoluteParticle::changeParticleState(int newState)
{

    KMCDebugger_MarkPre(particleStateName());


    KMCDebugger_Assert(m_totalParticles(particleState()), !=, 0, "trying to reduce particle type count below zero", info());

    m_totalParticles(particleState())--;

    m_totalParticles(newState)++;


    m_particleState = newState;

    KMCDebugger_PushImplication(this, particleStateName().c_str());

}


void SoluteParticle::reset()
{

    m_totalEnergy -= m_energy;

    setZeroEnergy();

    m_nNeighbors.zeros();

    m_nNeighborsSum = 0;

    for (Reaction * reaction : reactions())
    {
        reaction->reset();
    }


}

const uint &SoluteParticle::NX()
{
    return Site::NX();
}

const uint &SoluteParticle::NY()
{
    return Site::NY();
}

const uint &SoluteParticle::NZ()
{
    return Site::NZ();
}


KMCSolver *SoluteParticle::m_solver;

uvec4      SoluteParticle::m_totalParticles;

double     SoluteParticle::m_totalEnergy = 0;

#ifndef NDEBUG
particleSet SoluteParticle::m_affectedParticles = particleSet([] (SoluteParticle * s1, SoluteParticle * s2) {return s1->ID() < s2->ID();});
#else
particleSet SoluteParticle::m_affectedParticles;
#endif

uint SoluteParticle::ID_count = 0;
uint SoluteParticle::refCounter = 0;


ostream & operator << (ostream& os, const SoluteParticle& ss)
{
    os << ss.str();
    return os;
}
