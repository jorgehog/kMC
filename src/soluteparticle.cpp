#include "soluteparticle.h"

#include "kmcsolver.h"

#include "ignisinterface/solverevent.h"

#include "boundary/boundary.h"

#include "reactions/diffusion/tstdiffusion.h"

#include "reactions/diffusion/arrheniusdiffusion.h"

#include "potential/stressedsurface/stressedsurface.h"

#include "debugger/debugger.h"

#include <BADAss/badass.h>

using namespace kMC;


SoluteParticle::SoluteParticle(const uint species, bool sticky) :
    m_particleState(ParticleStates::solvant),
    m_site(NULL),
    m_x(UNSET_UINT),
    m_y(UNSET_UINT),
    m_z(UNSET_UINT),
    m_nNeighborsSum(0),
    m_energy(0),
    m_species(species),
    m_sticky(sticky),
    m_ID(ID_count++)
{

    BADAss(species, <, m_nSpecies, "invalid species.");

    if (!m_sticky)
    {
        initializeDiffusionReactions();
    }

    setVectorSizes();


    refCounter++;

}

SoluteParticle::~SoluteParticle()
{

    BADAss(refCounter, !=, 0);

    if (m_site != NULL)
    {
        disableSite();
    }

    clearAllReactions();

    m_nNeighbors.clear();


    refCounter--;

}


void SoluteParticle::setSite(const uint x, const uint y, const uint z)
{    

    BADAss(x, <, NX(), "mismatch in coordiantes. ", KMCBAI(info()));
    BADAss(y, <, NY(), "mismatch in coordiantes. ", KMCBAI(info()));
    BADAss(z, <, NZ(), "mismatch in coordiantes. ", KMCBAI(info()));

    trySite(x, y, z);

    BADAssBool(!m_site->isActive(), "particle already present at site.", KMCBAI( info()));

    m_site->associateWith(this);

    setupAllNeighbors();

    setNewParticleState(detectParticleState());

    forEachNeighborSiteDo_sendIndices([this] (Site *neighbor, uint i, uint j, uint k)
    {
        if (neighbor->isActive())
        {
            neighbor->associatedParticle()->addNeighbor(this, i, j, k);
        }

    });

    //tmp
    if (ss != NULL)
    {
        double dE = ss->evaluateFor(this);
        shiftEnergy(dE);
    }

    if (m_solver->solverEvent()->initialized())
    {

        forEachActiveReactionDo([] (Reaction *reaction)
        {
            reaction->forceUpdateFlag(Reaction::defaultUpdateFlag);
        });

        markAsAffected();

    }

    BADAss(m_site->associatedParticle(), ==, this, "mismatch in site and particle.");
    BADAssClose(getBruteForceTotalEnergy(), m_totalEnergy, 1E-3, "mismatch in total energy calculation.", [&] ()
    {
        cout << getBruteForceEnergy() << " " << m_totalEnergy << " " << getBruteForceEnergy() - m_totalEnergy << endl;
    });

    KMCDebugger_MarkPre("void");
    KMCDebugger_PushImplication(this, "enabled");
    KMCDebugger_MarkPartialStep("PARTICLE ACTIVATED");

}

void SoluteParticle::trySite(const uint x, const uint y, const uint z)
{
    m_x = x;
    m_y = y;
    m_z = z;

    BADAssBool(!Boundary::isBlocked(m_x, m_y, m_z), "coordinates were not set properly.");

    m_site = m_solver->getSite(x, y, z);
}

void SoluteParticle::resetSite()
{
    m_site = NULL;

    m_x = UNSET_UINT;
    m_y = UNSET_UINT;
    m_z = UNSET_UINT;

}

void SoluteParticle::disableSite()
{
    BADAss(m_site, !=, NULL, "disabeling disabled site. You need to call KMCSolver::despawnParticle to clean up particles.");
    BADAssBool(m_site->isActive(), "particle not present at site.", KMCBAI( info()));
    BADAss(m_site->associatedParticle(), ==, this, "mismatch in site and particle.");

    KMCDebugger_MarkPre(particleStateName());
    KMCDebugger_PushImplication(this, "disabled");

    m_site->desociate();
    m_site = NULL;

    clearNeighborhood();

    m_totalParticles(particleState())--;


    BADAssClose(getBruteForceTotalEnergy(), m_totalEnergy, 1E-3, "mismatch in total energy calculation.");

    KMCDebugger_MarkPartialStep("PARTICLE DISABLED");

}

void SoluteParticle::changePosition(const uint x, const uint y, const uint z)
{
    BADAssBool(x != m_x || y != m_y || z != m_z, "changing to the same site.");

    disableSite();

    setSite(x, y, z);

}

void SoluteParticle::setMainSolver(KMCSolver *solver)
{
    m_solver = solver;
}

void SoluteParticle::nSpecies(const uint _nSpecies, bool recalculatePotential)
{
    BADAss(_nSpecies, !=, 0);
    BADAssEqual(refCounter, 0);

    m_nSpecies = _nSpecies;

    if (recalculatePotential)
    {
        DiffusionReaction::setupPotential();
    }
}

void SoluteParticle::popAffectedParticle(SoluteParticle *particle)
{
    auto idx = m_affectedParticles.find(particle);

    if (idx == m_affectedParticles.end())
    {
        return;
    }

    BADAssBool(particle->isAffected());
    KMCDebugger_PopAffected(particle);

    m_affectedParticles.erase(idx);
}

void SoluteParticle::updateAffectedParticles()
{

    KMCDebugger_PushTraces();

    if (!solver()->localUpdating())
    {
        BADAssBool(m_affectedParticles.empty(), "non local updates should not use affected particles.", [] ()
        {
            for (SoluteParticle* particle : m_affectedParticles)
            {
                cout << *particle << endl;
            }
        });

        for (SoluteParticle *particle : solver()->particles())
        {
            if (particle != NULL)
            {
                particle->updateReactions();
            }
        }

    }

    else
    {

        for (SoluteParticle* particle : m_affectedParticles)
        {
            if (particle != NULL)
            {
                particle->updateReactions();
            }

        }

        clearAffectedParticles();

    }
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
bool SoluteParticle::isLegalToSpawn() const
{

    for (uint i = 0; i < 3; ++i)
    {
        for (uint j = 0; j < 3; ++j)
        {
            for (uint k = 0; k < 3; ++k)
            {
                if (i == j && j == k && k == 1)
                {
                    continue;
                }

                else if (diffusionReactions_fromIndex(i, j, k)->destinationSite() == NULL)
                {
                    continue;
                }

                else if (diffusionReactions_fromIndex(i, j, k)->destinationSite()->isActive())
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
    m_totalEnergy -= m_energy;

    m_energy = 0;

}

void SoluteParticle::setVectorSizes()
{

    m_nNeighbors.resize(Site::nNeighborsLimit());

}

void SoluteParticle::clearNeighborhood()
{
    forEachNeighborSiteDo_sendIndices([this] (Site *neighbor, uint i, uint j, uint k)
    {
        if (neighbor->isActive())
        {
            neighbor->associatedParticle()->removeNeighbor(this, i, j, k);
        }
    });

    setZeroEnergy();
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
    if (!solver()->localUpdating())
    {
        return;
    }

    m_affectedParticles.insert(this);
}

void SoluteParticle::updateReactions()
{

    for (Reaction* reaction : m_reactions)
    {
        if (reaction->isAllowed())
        {
            reaction->setRate();
            reaction->setLastUsedEnergy();
        }

        else
        {
            reaction->disable();
        }

        reaction->resetUpdateFlag();
    }
}

void SoluteParticle::setupAllNeighbors()
{

    BADAss(m_energy, ==, 0, "Particle energy should be cleared prior to this call.", KMCBAI( info()));

    m_nNeighbors.zeros();

    m_nNeighborsSum = 0;


    double dE = 0;

    forEachNeighborSiteDo_sendIndices([&dE, this] (Site * neighbor, uint i, uint j, uint k)
    {

        if (!neighbor->isActive())
        {
            return;
        }

        m_nNeighbors(Site::levelMatrix(i, j, k))++;

        m_nNeighborsSum++;


        dE += potentialBetween(neighbor->associatedParticle(), i, j, k);

    });

    shiftEnergy(dE);

}

void SoluteParticle::removeNeighbor(SoluteParticle *neighbor,
                                    const uint i,
                                    const uint j,
                                    const uint k)
{
    _updateNeighborProps(-1, neighbor, i, j, k);
}

void SoluteParticle::addNeighbor(SoluteParticle *neighbor,
                                 const uint i,
                                 const uint j,
                                 const uint k)
{
    _updateNeighborProps(+1, neighbor, i, j, k);
}

void SoluteParticle::_updateNeighborProps(const int sign,
                                          const SoluteParticle *neighbor,
                                          const uint i,
                                          const uint j,
                                          const uint k)
{
    KMCDebugger_MarkPre(neighbor->particleStateName());

    const uint &level = Site::levelMatrix(i, j, k);

    m_nNeighbors(level) += sign;

    m_nNeighborsSum += sign;

    double dE = sign*potentialBetween(neighbor, i, j, k);


    if (level == 0)
    {
        changeParticleState(detectParticleState());

        //tmp
        if (ss != NULL)
        {
            dE += ss->onNeighborChange(this, neighbor, i, j, k, sign);
        }

    }

    shiftEnergy(dE);

    if (m_solver->solverEvent()->initialized())
    {
        forEachActiveReactionDo([&level, &neighbor] (Reaction *reaction)
        {
            reaction->setDirectUpdateFlags(neighbor, level);
        });

        markAsAffected();

    }

#ifdef KMC_EXTRA_DEBUG
    BADAssClose(getBruteForceEnergy(), m_energy, 1E-3, "mismatch in energy calculation.");
#endif

    KMCDebugger_PushImplication(neighbor, neighbor->particleStateName().c_str());
}

void SoluteParticle::informOuterNeighbors() const
{
    SoluteParticle *shellOccupant;

    forShellDo(Site::nNeighborsLimit() + 1, [&shellOccupant, this] (Site *shellSite, int dx, int dy, int dz)
    {
        (void) dx;
        (void) dy;
        (void) dz;

        if (!shellSite->isActive())
        {
            return;
        }

        shellOccupant = shellSite->associatedParticle();

        shellOccupant->markAsAffected();
        shellOccupant->forEachActiveReactionDo([this] (Reaction *reaction)
        {
            reaction->setDirectUpdateFlags(this, Site::nNeighborsLimit());
        });
    });
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

    if (solver()->diffusionType() == KMCSolver::DiffusionTypes::None)
    {
        return;
    }

    BADAss(m_reactions.size(), ==, 0, "Sitereactions are already set", KMCBAI( info()));

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


                DiffusionReaction* diffusionReaction;


                if (solver()->diffusionType() == KMCSolver::DiffusionTypes::TST)
                {
                    diffusionReaction = new TSTDiffusion(this, dx, dy, dz);
                }

                else if (solver()->diffusionType() == KMCSolver::DiffusionTypes::Arrhenius)
                {
                    diffusionReaction = new ArrheniusDiffusion(this, dx, dy, dz);
                }

                else
                {
                    solver()->exit("Unknown diffusion type.");
                }


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

    else if (qualifiesAsSolvant())
    {
        return ParticleStates::solvant;
    }

    BADAssBool(qualifiesAsSurface());

    return ParticleStates::surface;

}


double SoluteParticle::potentialBetween(const SoluteParticle *other) const
{
    int X, Y, Z;

    distanceTo(other, X, Y, Z, true);

    X += Site::nNeighborsLimit();
    Y += Site::nNeighborsLimit();
    Z += Site::nNeighborsLimit();

    return potentialBetween(other, X, Y, Z);
}

double SoluteParticle::potentialBetween(const SoluteParticle *other,
                                        const uint i,
                                        const uint j,
                                        const uint k) const
{
    return DiffusionReaction::potential(i, j, k, m_species, other->species());
}

uint SoluteParticle::maxDistanceTo(const SoluteParticle *other) const
{
    return Site::maxDistanceBetween(m_x, m_y, m_z, other->x(), other->y(), other->z());
}




void SoluteParticle::setZeroTotalEnergy()
{
    m_totalEnergy = 0;

    BADAss(Site::solver()->particles().size(), ==, 0);
}

void SoluteParticle::setSticky(const bool sticky)
{
    if (m_sticky && !sticky)
    {
        initializeDiffusionReactions();
    }
    else if (sticky && !m_sticky)
    {
        clearAllReactions();
    }
    else
    {
        return;
    }

    updateReactions();

    solver()->reshuffleReactions();

    m_sticky = sticky;
}

uint SoluteParticle::nNeighborsSum() const
{
    BADAss(sum(m_nNeighbors), ==, m_nNeighborsSum, "Should be identical.", KMCBAI( info()));
    return m_nNeighborsSum;
}

double SoluteParticle::getBruteForceEnergy() const
{
    double energy = 0;

    forEachNeighborSiteDo([this, &energy] (Site *neighbor)
    {
        if (neighbor->isActive())
        {
            uint s = neighbor->associatedParticle()->species();
            int dx, dy, dz;
            distanceTo(neighbor->associatedParticle(),dx, dy, dz);
            double r = sqrt(dx*dx + dy*dy + dz*dz);
            energy += DiffusionReaction::strength(m_species, s)/pow(r, DiffusionReaction::rPower(m_species, s));

        }

    });

    //tmp
    if (ss != NULL)
    {
        double dE = ss->evaluateFor(this);
        energy += dE;
    }

    return energy;
}

void SoluteParticle::clearAll()
{
    BADAss(refCounter, ==, 0, "cannot clear static members with object alive.");

    //tmp
    if (ss != NULL)
    {
        delete ss;
        ss = NULL;
    }


    clearAffectedParticles();
    ID_count = 0;

    m_totalEnergy = 0;
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


    BADAss(m_totalParticles(particleState()), !=, 0, "trying to reduce particle type count below zero", KMCBAI( info()));

    m_totalParticles(particleState())--;

    m_totalParticles(newState)++;


    m_particleState = newState;

    KMCDebugger_PushImplication(this, particleStateName().c_str());

}

void SoluteParticle::shiftEnergy(const double amount)
{
    m_energy += amount;

    m_totalEnergy += amount;

}

double SoluteParticle::getBruteForceTotalEnergy()
{
    double totalEnergy = 0;

    for (const SoluteParticle *particle : m_solver->particles())
    {
        if (particle != NULL)
        {
            totalEnergy += particle->energy();
        }
    }

    return totalEnergy;

}


void SoluteParticle::reset()
{

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


Potential *SoluteParticle::ss = NULL;

KMCSolver *SoluteParticle::m_solver;

uvec4      SoluteParticle::m_totalParticles;

double     SoluteParticle::m_totalEnergy = 0;

function<bool(SoluteParticle *, SoluteParticle *)> SoluteParticle::compareFunc = [] (SoluteParticle * s1, SoluteParticle * s2) {return s1->ID() < s2->ID();};

#ifndef NDEBUG
particleSet SoluteParticle::m_affectedParticles = particleSet(SoluteParticle::compareFunc);
#else
particleSet SoluteParticle::m_affectedParticles;
#endif

uint        SoluteParticle::m_nSpecies = 1;


uint SoluteParticle::ID_count = 0;
uint SoluteParticle::refCounter = 0;

ostream & operator << (ostream& os, const SoluteParticle& ss)
{
    os << ss.str();
    return os;
}
