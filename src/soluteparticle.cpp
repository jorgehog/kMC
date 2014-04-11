#include "soluteparticle.h"

#include "kmcsolver.h"

#include "debugger/debugger.h"

#include "reactions/diffusion/diffusionreaction.h"


using namespace kMC;


SoluteParticle::SoluteParticle() :
    m_particleState(ParticleStates::solvant),
    m_nNeighborsSum(0),
    m_energy(0),
    m_ID(refCounter)
{

    initializeDiffusionReactions();

    setVectorSizes();


    refCounter++;

}

SoluteParticle::~SoluteParticle()
{

    KMCDebugger_Assert(refCounter, !=, 0);

    removeCurrentSite();

    clearAllReactions();


    refCounter--;

}


void SoluteParticle::setSite(kMC::Site *site)
{

    m_site = site;


    markAsAffected();

    informNeighborhoodOnAddition();

    cout << "Derp: have to update neighbor list here" << endl;
    m_site->activate();

    //DERP: do this checks at once?
    if (qualifiesAsCrystal())
    {
        setNewParticleState(ParticleStates::crystal);
    }

    else if (qualifiesAsSurface())
    {
        setNewParticleState(ParticleStates::surface);
    }

    else
    {
        setNewParticleState(ParticleStates::solvant);
    }

    setNeighboringDirectUpdateFlags();

    for (Reaction * reaction : m_reactions)
    {
        reaction->setDirectUpdateFlags(this);
    }


    KMCDebugger_MarkPartialStep("PARTICLE SITE CHANGE END");

}

void SoluteParticle::trySite(Site *site)
{
    m_site = site;
}

void SoluteParticle::removeCurrentSite()
{

    setNeighboringDirectUpdateFlags();

    informNeighborhoodOnRemoval();

    m_totalParticles(particleState())--;

    m_totalEnergy -= m_energy;

    m_site->deactivate();

}

void SoluteParticle::changeSite(Site *newSite)
{
    removeCurrentSite();

    setSite(newSite);

}

void SoluteParticle::loadConfig(const Setting &setting)
{
    setInitialNNeighborsToCrystallize(getSurfaceSetting<uint>(setting, "nNeighboursToCrystallize"));
}

void SoluteParticle::popAffectedParticle(SoluteParticle *particle)
{
    KMCDebugger_AssertBool(particle->isAffected());
    KMCDebugger_PopAffected(particle);

    m_affectedParticles.erase(m_affectedParticles.find(particle));
}

void SoluteParticle::updateAffectedParticles()
{

    KMCDebugger_PushTraces();

    for (SoluteParticle* particle : m_affectedParticles)
    {
        particle->updateReactions();
    }

    clearAffectedParticles();

}





void SoluteParticle::setParticleState(int newState)
{

    m_totalParticles(m_particleState)--;
    m_totalParticles(newState)++;

    m_particleState = newState;

}




bool SoluteParticle::qualifiesAsCrystal()
{

    if (nNeighbors() == 26)
    {
        return true;
    }

    return false;

}

bool SoluteParticle::qualifiesAsSurface()
{
    return m_site->hasNeighboring(ParticleStates::crystal) && !qualifiesAsCrystal();
}

//All reactions must be legal if site is allowed to spawn.
bool SoluteParticle::isLegalToSpawn()
{

    for (Reaction * r : m_reactions)
    {
        if (!r->isAllowed())
        {
            return false;
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

const string SoluteParticle::info(int xr, int yr, int zr, string desc) const
{
    stringstream s;

    s << "SoluteParticle@";
    s << m_site->info(xr, yr, zr, desc);

    return s.str();

}

void SoluteParticle::setZeroEnergy()
{
    KMCDebugger_Assert(m_nNeighborsSum, ==, 0, "Energy is not zero.", info());
    m_energy = 0;

}

void SoluteParticle::setVectorSizes()
{

    m_neighboringParticles.resize(Site::nNeighborsLimit());

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

    m_neighboringParticles.resize(Site::nNeighborsLimit());

    m_nNeighbors.zeros();

    m_nNeighborsSum = 0;

    m_totalEnergy -= m_energy;

    m_energy = 0;


    for (vector<SoluteParticle*> & neighborShell : m_neighboringParticles)
    {
        neighborShell.clear();
    }


    double dE;

    uint level;

    site()->forEachNeighborDo_sendIndices([&level, this] (Site * neighbor, uint i, uint j, uint k)
    {

        if (!neighbor->isActive())
        {
            return;
        }

        level = Site::levelMatrix(i, j, k);

        m_neighboringParticles.at(level).push_back(neighbor->getAssociatedParticle());

        m_nNeighbors(level)++;

        m_nNeighborsSum++;

        dE = DiffusionReaction::potential(i,  j,  k);

        m_energy += dE;

        m_totalEnergy += dE;

    });



}

void SoluteParticle::removeNeighbor(SoluteParticle *neighbor, uint level)
{
    m_neighboringParticles.at(level).erase(std::find(m_neighboringParticles.at(level).begin(),
                                                     m_neighboringParticles.at(level).end(),
                                                     neighbor),
                                           m_neighboringParticles.at(level).end());
}

void SoluteParticle::addNeighbor(SoluteParticle *neighbor, uint level)
{
    m_neighboringParticles.at(level).push_back(neighbor);
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


double SoluteParticle::potentialBetween(const SoluteParticle *other)
{
    int X, Y, Z;

    m_site->distanceTo(other->site(), X, Y, Z, true);

    X += Site::nNeighborsLimit();
    Y += Site::nNeighborsLimit();
    Z += Site::nNeighborsLimit();

    return DiffusionReaction::potential(X, Y, Z);
}

void SoluteParticle::setNeighboringDirectUpdateFlags()
{

    for (SoluteParticle *neighbor : m_neighboringParticles)
    {
        {
            //This approach assumes that recursive updating of non-neighboring sites
            //WILL NOT ACTIVATE OR DEACTIVATE any sites, simply change their state,
            //and thus not interfere with any flags set here, not require flags of their own.
            for (Reaction * reaction : neighbor->reactions())
            {
                reaction->setDirectUpdateFlags(this);
            }

            neighbor->markAsAffected();

        }
    }

}

void SoluteParticle::setZeroTotalEnergy()
{
    m_totalEnergy = 0;

    KMCDebugger_Assert(solver()->particles().size(), ==, 0);
}

void SoluteParticle::clearAll()
{

    clearAffectedParticles();

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

    return static_cast<double>(nSolutionParticles())/(NX()*NY()*NZ() - nCrystals());
}

double SoluteParticle::getCurrentRelativeCrystalOccupancy()
{
    return static_cast<double>(nCrystals())/(NX()*NY()*NZ());
}

void SoluteParticle::setInitialNNeighborsToCrystallize(const uint &nNeighborsToCrystallize)
{
    m_nNeighborsToCrystallize = nNeighborsToCrystallize;
}


void SoluteParticle::resetNNeighborsToCrystallizeTo(const uint &nNeighborsToCrystallize)
{
    setInitialNNeighborsToCrystallize(nNeighborsToCrystallize);
}

void SoluteParticle::setNewParticleState(int newState)
{

    KMCDebugger_MarkPre(particleStateName());


    KMCDebugger_Assert(m_totalParticles(particleState()), !=, 0, "trying to reduce particle type count below zero", info());

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
    //    m_totalParticles(particleState())--;
    //        m_totalDeactiveParticles(particleState())++;

    //    if (isFixedCrystalSeed())
    //    {
    //        m_isFixedCrystalSeed = false;

    ////        initializeDiffusionReactions();
    //    }

    //    m_cannotCrystallize = false;

    //    m_totalDeactiveParticles(particleState())--;

    //    m_particleState = ParticleStates::solvant;

    //    m_totalDeactiveParticles(ParticleStates::solution)++;


    m_totalEnergy -= m_energy;

    setZeroEnergy();

    m_nNeighbors.zeros();

    m_nNeighborsSum = 0;

    for (Reaction * reaction : reactions())
    {
        reaction->reset();
    }


}

void SoluteParticle::queueAffectedParticles()
{

    for (SoluteParticle *neighbor : m_neighboringParticles)
    {

        for (Reaction * reaction : neighbor->reactions())
        {
            reaction->registerUpdateFlag(Reaction::defaultUpdateFlag);
        }

        neighbor->markAsAffected();
    }

}



void SoluteParticle::informNeighborhoodOnAddition()
{

    uint level;
    double dE;

    forEachNeighborDo_sendIndices([&dE, &level, this] (Site *neighbor, uint i, uint j, uint k)
    {

        level = m_levelMatrix(i, j, k);


        KMCDebugger_AssertBool(!((change < 0) && (neighbor->nNeighborsSum() == 0)), "Call initiated to set negative nNeighbors.", neighbor->info());
        KMCDebugger_AssertBool(!((change < 0) && (neighbor->nNeighbors(level) == 0)), "Call initiated to set negative neighbor.", neighbor->info());


        neighbor->m_nNeighbors(level)+= change;
        neighbor->m_nNeighborsSum += change;

        dE = change*DiffusionReaction::potential(i,  j,  k);

        neighbor->m_energy += dE;

        m_totalEnergy += dE;

    });

    for (uint level = 0; level < Site::nNeighborsLimit(); ++level)
    {
        for (SoluteParticle *neighbor : m_neighbors.at(level))
        {
            neighbor->removeNeighbor(this);
        }
    }


}


void SoluteParticle::informNeighborhoodOnRemoval()
{

    uint level;
    double dE;

    forEachNeighborDo_sendIndices([&change, &dE, &level, this] (Site *neighbor, uint i, uint j, uint k)
    {

        level = m_levelMatrix(i, j, k);


        KMCDebugger_AssertBool(!((change < 0) && (neighbor->nNeighborsSum() == 0)), "Call initiated to set negative nNeighbors.", neighbor->info());
        KMCDebugger_AssertBool(!((change < 0) && (neighbor->nNeighbors(level) == 0)), "Call initiated to set negative neighbor.", neighbor->info());


        neighbor->m_nNeighbors(level)+= change;
        neighbor->m_nNeighborsSum += change;


        dE = change*DiffusionReaction::potential(i,  j,  k);

        neighbor->m_energy += dE;

        m_totalEnergy += dE;

    });


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



uvec4      SoluteParticle::m_totalParticles;

double     SoluteParticle::m_totalEnergy = 0;


//Don't panic: Just states that the sets pointers should be sorted by the site objects and not their random valued pointers.
set<SoluteParticle*, function<bool(SoluteParticle*, SoluteParticle*)> > SoluteParticle::m_affectedParticles = set<SoluteParticle*, function<bool(SoluteParticle*, SoluteParticle*)> >([] (SoluteParticle * s1, SoluteParticle * s2) {return s1->ID() < s2->ID();});


uint SoluteParticle::refCounter = 0;


ostream & operator << (ostream& os, const SoluteParticle& ss)
{
    os << ss.str();
    return os;
}
