#include "soluteparticle.h"

#include "kmcsolver.h"

#include "debugger/debugger.h"

#include "reactions/diffusion/diffusionreaction.h"


using namespace kMC;


SoluteParticle::SoluteParticle() :
    m_particleState(ParticleStates::solvant),
    m_ID(refCounter)
{

    initializeDiffusionReactions();

    refCounter++;
}

SoluteParticle::~SoluteParticle()
{
    m_totalParticles(particleState())--;
    //    m_totalDeactiveParticles(particleState())++;

    clearAllReactions();

    m_site->deactivate();

    KMCDebugger_Assert(refCounter, !=, 0);

    refCounter--;

}


void SoluteParticle::setSite(kMC::Site *site)
{

    m_site = site;


    markAsAffected();

    site->informNeighborhoodOnChange(+1);

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

    m_site->informNeighborhoodOnChange(-1);

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

    if (m_site->nNeighbors() >= m_nNeighborsToCrystallize)
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


//void SoluteParticle::activate()
//{

//    flipActive();


//    if (isSurface())
//    {
//        setParticleState(ParticleStates::crystal);
//    }
//    else
//    {
//        KMCDebugger_PushImplication(this, "activation");
//    }


//    setNeighboringDirectUpdateFlags();

//    for (Reaction * reaction : m_reactions)
//    {
//        reaction->setDirectUpdateFlags(this);
//    }


//    KMCDebugger_MarkPartialStep("ACTIVATION COMPLETE");

//    m_totalActiveSites++;

//    KMCDebugger_AssertBool(!(isSurface() && isActive()), "surface should not be active.", info());

//}


//void SoluteParticle::deactivate()
//{

//    if (isFixedCrystalSeed())
//    {
//        deactivateFixedCrystal();
//        return;
//    }

//    flipDeactive();

//    //if we deactivate a crystal site, we have to potentially
//    //reduce the surface by adding more sites as solution sites.
//    //Site will change only if it is not surrounded by any crystals.
//    if (isCrystal())
//    {
//        KMCDebugger_Assert(particleState(), !=, ParticleStates::surface);
//        KMCDebugger_Assert(particleState(), !=, ParticleStates::solvant);

//        setParticleState(ParticleStates::surface);
//    }


//    else if (qualifiesAsSurface())
//    {

//        KMCDebugger_Assert(particleState(), ==, ParticleStates::solvant);

//        setNewParticleState(ParticleStates::surface);

//    }

//    else
//    {
//        KMCDebugger_PushImplication(this, "deactivation");
//    }


//    setNeighboringDirectUpdateFlags();

//    KMCDebugger_MarkPartialStep("DEACTIVATION COMPLETE");

//    m_totalActiveSites--;

//}

void SoluteParticle::clearAll()
{
    m_nNeighborsToCrystallize = KMCSolver::UNSET_UINT;

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

uint       SoluteParticle::m_nNeighborsToCrystallize = KMCSolver::UNSET_UINT;


uvec4      SoluteParticle::m_totalParticles;


//Don't panic: Just states that the sets pointers should be sorted by the site objects and not their random valued pointers.
set<SoluteParticle*, function<bool(SoluteParticle*, SoluteParticle*)> > SoluteParticle::m_affectedParticles = set<SoluteParticle*, function<bool(SoluteParticle*, SoluteParticle*)> >([] (SoluteParticle * s1, SoluteParticle * s2) {return s1->ID() < s2->ID();});


uint SoluteParticle::refCounter = 0;


ostream & operator << (ostream& os, const SoluteParticle& ss)
{
    os << ss.str();
    return os;
}
