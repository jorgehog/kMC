#include "kmcsolver.h"

#include "soluteparticle.h"

#include "reactions/reaction.h"
#include "reactions/diffusion/diffusionreaction.h"

#include "boundary/boundary.h"

#include "ignisinterface/solverevent.h"
#include "ignisinterface/kmcparticles.h"

#include <sys/time.h>

#include <armadillo>

#include <iostream>
#include <iomanip>
#include <fstream>


using namespace arma;
using namespace std;
using namespace kMC;

KMCSolver::KMCSolver(const Setting & root)
{

    onConstruct();

    const Setting & SystemSettings = getSurfaceSetting(root, "System");
    const Setting & SolverSettings = getSurfaceSetting(root, "Solver");
    const Setting & diffusionSettings = getSetting(root, {"Reactions", "Diffusion"});


    Reaction::loadConfig(getSurfaceSetting(root, "Reactions"));

    DiffusionReaction::loadConfig(diffusionSettings);

    Site::loadConfig(SystemSettings);

    SoluteParticle::loadConfig(SystemSettings);


    setNumberOfCycles(
                getSurfaceSetting<uint>(SolverSettings, "nCycles"));

    setCyclesPerOutput(
                getSurfaceSetting<uint>(SolverSettings, "cyclesPerOutput"));

    setRNGSeed(
                getSurfaceSetting<uint>(SolverSettings, "seedType"),
                getSurfaceSetting<int>(SolverSettings, "specificSeed"));

    setTargetConcentration(
                getSurfaceSetting<double>(SystemSettings, "SaturationLevel"));


    uvec3 boxSize;

    boxSize(0) = getSurfaceSetting(SystemSettings, "BoxSize")[0];
    boxSize(1) = getSurfaceSetting(SystemSettings, "BoxSize")[1];
    boxSize(2) = getSurfaceSetting(SystemSettings, "BoxSize")[2];

    setBoxSize(boxSize);

}

KMCSolver::KMCSolver()
{
    onConstruct();
}

KMCSolver::~KMCSolver()
{

    finalizeObject();

    clearSites();

    Site::clearAll();
    DiffusionReaction::clearAll();
    Boundary::clearAll();

    refCounter--;

}

void KMCSolver::checkRefCounter()
{
    if (refCounter > 1)
    {
        cerr << "ERROR: "<< refCounter << " solver objects alive when freeing.\n";
        cerr << "Static member variables of objects IN USE by living solver WILL BE FREED." << endl;

        exit();
    }
}

void KMCSolver::onConstruct()
{

    m_NX = UNSET_UINT;
    m_NY = UNSET_UINT;
    m_NZ = UNSET_UINT;

    m_targetConcentration = 0;

    totalTime = 0;

    cycle = 1;

    m_nCycles = 0;

    outputCounter = 0;

    m_kTot = 0;


    Boundary::setMainSolver(this);

    Reaction::setMainSolver(this);

    Site::setMainSolver(this);

    KMCParticles *particles = new KMCParticles(this);
    MainLattice::setCurrentParticles(*particles);

    m_mainLattice = new MainLattice();

    SolverEvent *solverEvent = new SolverEvent(this);
    m_mainLattice->addEvent(*solverEvent);

    refCounter++;

}

void KMCSolver::finalizeObject()
{
    checkRefCounter();

    for (SoluteParticle *particle : m_particles)
    {
        delete particle;
    }

    KMCDebugger_Finalize();

    //tmp
    delete m_mainLattice;
    Event<uint>::resetEventParameters();

}

void KMCSolver::mainloop()
{
    m_mainLattice->eventLoop(m_nCycles);
}


void KMCSolver::reset()
{

    finalizeObject();


    totalTime = 0;

    cycle = 1;

    outputCounter = 0;

    m_kTot = 0;


    m_allPossibleReactions.clear();

    m_accuAllRates.clear();

    m_availableReactionSlots.clear();


    SoluteParticle::clearAffectedParticles();

    forEachSiteDo([] (Site * site)
    {
        site->reset();
    });

    KMCDebugger_Assert(accu(SoluteParticle::totalParticlesVector()), ==, 0);

    //    KMCDebugger_Assert(Site::totalDeactiveParticles(ParticleStates::solution), ==, m_NX*m_NY*m_NZ);

    KMCDebugger_AssertClose(SoluteParticle::totalEnergy(), 0, 1E-5);

    KMCDebugger_Assert(SoluteParticle::nParticles(), ==, 0);

    SoluteParticle::setZeroTotalEnergy();

    Site::finalizeBoundaries();

    setRNGSeed(Seed::specific, Seed::initialSeed);

    Site::initializeBoundaries();

    KMCDebugger_Init();


    //TMP
    KMCParticles *particles = new KMCParticles(this);
    MainLattice::setCurrentParticles(*particles);

    m_mainLattice = new MainLattice();

    SolverEvent *solverEvent = new SolverEvent(this);
    m_mainLattice->addEvent(*solverEvent);

}



void KMCSolver::dumpXYZ()
{
    cout << "Storing XYZ: " << outputCounter << endl;

    stringstream s;
    s << "kMC" << outputCounter++ << ".xyz";

    ofstream o;
    o.open("outfiles/" + s.str());



    stringstream surface;
    stringstream crystal;
    stringstream solution;

    uint nLines = 0;
    s.str(string());

    for (SoluteParticle *particle : m_particles)
    {

        s << "\n"
          << particle->particleStateShortName() << " "
          << particle->x() << " " << particle->y() << " " << particle->z() << " "
          << particle->nNeighborsSum() << " "
          << particle->energy();

        if (particle->isSurface())
        {
            surface << s.str();
        }

        else if (particle->isCrystal())
        {
            crystal << s.str();
        }

        else
        {
            solution << s.str();
        }

        s.str(string());
        nLines++;

    }

    o << nLines << "\n - " << surface.str() << crystal.str() << solution.str();
    o.close();

}

void KMCSolver::onAllRatesChanged()
{
    //Change if added reactions are not on form ..exp(beta...)

    m_kTot = 0;


    uint i = 0;
    for (Reaction * r : m_allPossibleReactions)
    {

        KMCDebugger_Assert(r->rate(), !=, Reaction::UNSET_RATE, "Rates should not be all changed before they are set once.");
        KMCDebugger_Assert(r->rate(), >, 0, "Rate should be positive.");
        KMCDebugger_AssertBool(!r->hasVacantStatus(), "Reaction should be enabled.");

        m_kTot += r->rate();

        m_accuAllRates.at(r->address()) = m_kTot;

        i++;

    }

}

void KMCSolver::registerReactionChange(Reaction *reaction, const double &newRate)
{

    const double & prevRate = reaction->rate();

    if (prevRate == newRate)
    {
        return;
    }

    else if (prevRate == Reaction::UNSET_RATE)
    {

        KMCDebugger_Assert(newRate, !=, Reaction::UNSET_RATE);

        m_kTot += newRate;


        //If there is a vacancy, we simply fill it.
        if (!m_availableReactionSlots.empty())
        {

            const uint slot = m_availableReactionSlots.at(m_availableReactionSlots.size() - 1);

            reaction->setAddress(slot);

            m_allPossibleReactions.at(slot) = reaction;

            m_accuAllRates.at(slot) = prevAccuAllRatesValue(slot);


            updateAccuAllRateElements(slot, m_accuAllRates.size(), newRate);


            m_availableReactionSlots.pop_back();

        }

        //if not, we make a new element
        else
        {
            m_allPossibleReactions.push_back(reaction);

            m_accuAllRates.push_back(m_kTot);

            reaction->setAddress(m_allPossibleReactions.size()-1);
        }

    }

    else if (newRate == Reaction::UNSET_RATE)
    {

        KMCDebugger_AssertBool(!reaction->isAllowed(), "Allowed reaction set to unset rate.");

        m_kTot -= prevRate;


        KMCDebugger_Assert(reaction->address(), !=, Reaction::UNSET_ADDRESS);
        KMCDebugger_AssertBool(!isEmptyAddress(reaction->address()), "address is already set as empty.");

        m_availableReactionSlots.push_back(reaction->address());

        //reset the accuallrates value to the previous value or zero, so that when we swap in a reaction to a vacant spot,
        //we simply add the rate value on top of all higher elements.
        m_accuAllRates.at(reaction->address()) = prevAccuAllRatesValue(reaction->address());

        updateAccuAllRateElements(reaction->address() + 1, m_accuAllRates.size(), -prevRate);

        reaction->setAddress(Reaction::UNSET_ADDRESS);

    }

    else
    {
        double deltaRate = (newRate - prevRate);

        m_kTot += deltaRate;


        updateAccuAllRateElements(reaction->address(), m_accuAllRates.size(), deltaRate);

    }

}

void KMCSolver::reshuffleReactions()
{

    uint nVacancies = m_availableReactionSlots.size();

    uint firstVacancy;
    uint lastReaction = m_allPossibleReactions.size() - 1;

    uint numberOfSwaps = 0;
    uint trailingVacancies   = 0;

    std::sort(m_availableReactionSlots.begin(), m_availableReactionSlots.end());

    //While we have not yet filled all vacancies
    while (trailingVacancies < nVacancies)
    {

        firstVacancy = m_availableReactionSlots.at(numberOfSwaps);

        //(trailingVacancies - numberOfSwaps) is the number of additional shifts we need to make away
        //from the last vacant spot. This is greater than zero only if we have trailing vacant sites
        //present before swapping.
        while (lastReaction == m_availableReactionSlots.at((nVacancies - 1) - (trailingVacancies - numberOfSwaps)))
        {
            lastReaction--;
            trailingVacancies++;

            //This terminates the function.
            if (trailingVacancies == nVacancies)
            {
                postReactionShuffleCleanup(nVacancies);

                return;
            }
        }


        swapReactionAddresses(firstVacancy, lastReaction);

        lastReaction--;
        trailingVacancies++;

        numberOfSwaps++;

    }

    postReactionShuffleCleanup(nVacancies);

}

void KMCSolver::swapReactionAddresses(const uint dest, const uint orig)
{

    KMCDebugger_AssertBool(isEmptyAddress(dest),  "destination should be empty.");
    KMCDebugger_AssertBool(!isEmptyAddress(orig), "origin should not be empty.");

    Reaction * swappedReaction = m_allPossibleReactions.at(orig);
    Reaction * oldReaction     = m_allPossibleReactions.at(dest);

    KMCDebugger_Assert(orig,                    ==, swappedReaction->address(), "mismatch in address.", swappedReaction->getFinalizingDebugMessage());
    KMCDebugger_Assert(Reaction::UNSET_ADDRESS, ==, oldReaction->address(),     "mismatch in address.", oldReaction->getFinalizingDebugMessage());

    KMCDebugger_AssertBool(swappedReaction->isAllowed(), "swapped reaction should be allowed and active.", swappedReaction->getFinalizingDebugMessage());
    KMCDebugger_AssertBool(!oldReaction->isAllowed()   , "old reaction should not be allowed.",            oldReaction->getFinalizingDebugMessage());

    m_allPossibleReactions.at(dest) = swappedReaction;

    swappedReaction->setAddress(dest);
    oldReaction->setAddress(Reaction::UNSET_ADDRESS);

    updateAccuAllRateElements(dest, orig, swappedReaction->rate());



}

void KMCSolver::postReactionShuffleCleanup(const uint nVacancies)
{

    //Optimize further: Do not use resize, but rather keep the limit in memory.
    m_allPossibleReactions.resize(m_allPossibleReactions.size() - nVacancies);
    m_accuAllRates.resize(m_accuAllRates.size() - nVacancies);

    m_availableReactionSlots.clear();

    KMCDebugger_Assert(m_allPossibleReactions.size(), ==, m_accuAllRates.size(), "These vectors should be equal of length.");
    m_accuAllRates.size() == 0
            ? (void) 0
            : KMCDebugger_AssertClose(m_accuAllRates.at(m_accuAllRates.size()- 1),
                                      m_kTot, minRateThreshold(),
                                      "kTot should be the last element of accuAllRates");

}


bool KMCSolver::isEmptyAddress(const uint address)
{
    return std::find(m_availableReactionSlots.begin(), m_availableReactionSlots.end(), address)
            != m_availableReactionSlots.end();
}

string KMCSolver::getReactionVectorDebugMessage()
{
    stringstream s;

    s << "vacant addresses: \n";

    for (uint addr : m_availableReactionSlots)
    {
        s << addr << "\n";
    }

    s << "\npossible reactions: \n";

    for (Reaction * r : m_allPossibleReactions)
    {
        s << r->str() << " " << r->propertyString() << "\n";
    }

    return s.str();

}

void KMCSolver::dumpOutput()
{

    cout << setw(5) << right << setprecision(1) << fixed
         << (double)cycle/m_nCycles*100 << "%   "
         << outputCounter
         << "\n";
    cout << setprecision(6);
}


void KMCSolver::forEachSiteDo(function<void (Site *)> applyFunction) const
{
    for (uint x = 0; x < m_NX; ++x)
    {
        for (uint y = 0; y < m_NY; ++y)
        {
            for (uint z = 0; z < m_NZ; ++z)
            {
                applyFunction(sites[x][y][z]);
            }
        }
    }
}

void KMCSolver::forEachSiteDo_sendIndices(function<void (Site *, uint, uint, uint)> applyFunction) const
{
    for (uint x = 0; x < m_NX; ++x)
    {
        for (uint y = 0; y < m_NY; ++y)
        {
            for (uint z = 0; z < m_NZ; ++z)
            {
                applyFunction(sites[x][y][z], x, y, z);
            }
        }
    }
}

void KMCSolver::forEachActiveSiteDo(function<void (Site *)> applyFunction) const
{
    for (uint x = 0; x < m_NX; ++x)
    {
        for (uint y = 0; y < m_NY; ++y)
        {
            for (uint z = 0; z < m_NZ; ++z)
            {
                if (sites[x][y][z]->isActive())
                {
                    applyFunction(sites[x][y][z]);
                }
            }
        }
    }
}

void KMCSolver::forEachActiveSiteDo_sendIndices(function<void (Site *, uint, uint, uint)> applyFunction) const
{
    for (uint x = 0; x < m_NX; ++x)
    {
        for (uint y = 0; y < m_NY; ++y)
        {
            for (uint z = 0; z < m_NZ; ++z)
            {
                if (sites[x][y][z]->isActive())
                {
                    applyFunction(sites[x][y][z], x, y, z);
                }
            }
        }
    }
}


void KMCSolver::initializeSites()
{

    sites = new Site***[m_NX];

    for (uint x = 0; x < m_NX; ++x)
    {
        sites[x] = new Site**[m_NY];

        for (uint y = 0; y < m_NY; ++y)
        {
            sites[x][y] = new Site*[m_NZ];

            for (uint z = 0; z < m_NZ; ++z)
            {
                sites[x][y][z] = new Site(x, y, z);
            }
        }
    }


    initializeSiteNeighborhoods();

}


void KMCSolver::clearSites()
{

    KMCDebugger_SetEnabledTo(false);

    for (uint i = 0; i < m_NX; ++i)
    {
        for (uint j = 0; j < m_NY; ++j)
        {
            for (uint k = 0; k < m_NZ; ++k)
            {
                delete sites[i][j][k];
            }

            delete [] sites[i][j];
        }

        delete [] sites[i];
    }

    delete [] sites;


    KMCDebugger_Assert(accu(SoluteParticle::totalParticlesVector()), ==, 0);

    KMCDebugger_AssertClose(SoluteParticle::totalEnergy(), 0, 1E-5);

    SoluteParticle::clearAffectedParticles();
    SoluteParticle::setZeroTotalEnergy();

    Reaction::clearAll();

    m_allPossibleReactions.clear();

    m_accuAllRates.clear();

    m_availableReactionSlots.clear();


    KMCDebugger_ResetEnabled();

}

void KMCSolver::setBoxSize_KeepSites(const uvec3 &boxSizes)
{

    // have to clear sites and all that jazz because of vector setup...

}


bool KMCSolver::spawnParticle(SoluteParticle *particle, uint x, uint y, uint z, bool checkIfLegal)
{
    spawnParticle(particle, getSite(x, y, z), checkIfLegal);
}


bool KMCSolver::spawnParticle(SoluteParticle *particle, Site *site, bool checkIfLegal)
{

    KMCDebugger_AssertBool(!(site->isActive() && !checkIfLegal), "spawning particle on top of another");

    if (site->isActive())
    {
        return false;
    }

    particle->trySite(site);

    if (checkIfLegal)
    {
        if (!particle->isLegalToSpawn())
        {
            return false;
        }
    }

    particle->setSite(site);

    m_particles.push_back(particle);

    KMCDebugger_PushTraces();

    return true;

}

void KMCSolver::forceSpawnParticle(Site *site)
{
    SoluteParticle *particle = new SoluteParticle();

    spawnParticle(particle, site, false);
}

void KMCSolver::despawnParticle(Site *site)
{

    KMCDebugger_AssertBool(site->isActive());

    uint i = 0;

    for (SoluteParticle *particle : m_particles)
    {
        if (particle->site() == site)
        {
            break;
        }
        i++;
    }

    KMCDebugger_Assert(i, !=, m_particles.size());

    delete m_particles.at(i);

    m_particles.erase(m_particles.begin() + i);

}


void KMCSolver::initializeCrystal(const double relativeSeedSize)
{

    if (relativeSeedSize > 1.0)
    {
        cerr << "The seed size cannot exceed the box size." << endl;
        KMCSolver::exit();
    }

    else if (relativeSeedSize < 0)
    {
        cerr << "The seed size cannot be negative." << endl;
        KMCSolver::exit();
    }

    KMCDebugger_SetEnabledTo(false);

    uint crystalSizeX = round(m_NX*relativeSeedSize);
    uint crystalSizeY = round(m_NY*relativeSeedSize);
    uint crystalSizeZ = round(m_NZ*relativeSeedSize);

    uint crystalStartX = m_NX/2 - crystalSizeX/2;
    uint crystalStartY = m_NY/2 - crystalSizeY/2;
    uint crystalStartZ = m_NZ/2 - crystalSizeZ/2;

    uint crystalEndX = crystalStartX + crystalSizeX;
    uint crystalEndY = crystalStartY + crystalSizeY;
    uint crystalEndZ = crystalStartZ + crystalSizeZ;

    for (uint i = 0; i < m_NX; ++i)
    {
        for (uint j = 0; j < m_NY; ++j)
        {
            for (uint k = 0; k < m_NZ; ++k)
            {


                if (i >= crystalStartX && i < crystalEndX)
                {
                    if (j >= crystalStartY && j < crystalEndY)
                    {
                        if (k >= crystalStartZ && k < crystalEndZ)
                        {
                            if (!((i == m_NX/2 && j == m_NY/2 && k == m_NZ/2)))
                            {
                                SoluteParticle * particle = new SoluteParticle();
                                (void)spawnParticle(particle, i, j, k, false);
                            }

                        }
                    }
                }

            }
        }
    }

    initializeSolutionBath();

    KMCDebugger_ResetEnabled();

}

void KMCSolver::initializeSolutionBath()
{

    uint x, y, z, N, n;
    bool spawned;


    N = (NX()*NY()*NZ() - SoluteParticle::nParticles())*targetConcentration();

    n = 0;

    while (n != N)
    {


        SoluteParticle *particle = new SoluteParticle();

        spawned = false;

        while (!spawned)
        {

            x = KMC_RNG_UNIFORM()*NX();

            y = KMC_RNG_UNIFORM()*NY();

            z = KMC_RNG_UNIFORM()*NZ();

            spawned = spawnParticle(particle, x, y, z, true);

        }

        n++;
    }
}




void KMCSolver::getRateVariables()
{

    SoluteParticle::updateAffectedParticles();

    reshuffleReactions();

}


uint KMCSolver::getReactionChoice(double R)
{

    KMCDebugger_Assert(m_accuAllRates.size(), !=, 0, "No active reactions.");

    uint imax = m_accuAllRates.size() - 1;
    uint MAX = imax;
    uint imin = 0;
    uint imid = 1;

    // continue searching while imax != imin + 1
    while (imid != imin)
    {

        // calculate the midpoint for roughly equal partition
        imid = imin + (imax - imin)/2;

        //Is the upper limit above mid?
        if (R > m_accuAllRates.at(imid))
        {

            if (imid == MAX)
            {
                return MAX;
            }

            //Are we just infront of the limit?
            else if (R < m_accuAllRates.at(imid + 1))
            {
                //yes we were! Returning current mid + 1.
                //If item i in accuAllrates > R, then reaction i is selected.
                //This because there is no zero at the start of accuAllrates.

                return imid + 1;
            }

            //No we're not there yet, so we search above us.
            else
            {
                imin = imid + 1;
            }
        }

        //No it's not. Starting new search below mid!
        else
        {

            if (imid == 0)
            {
                return 0;
            }

            imax = imid;
        }


    }

    return imid + 1;

}

void KMCSolver::setBoxSize(const uvec3 boxSize, bool check, bool keepSystem)
{

    if (keepSystem)
    {
        setBoxSize_KeepSites(boxSize);
        return;
    }

    if (m_NX != UNSET_UINT && m_NY != UNSET_UINT && m_NZ != UNSET_UINT)
    {
        clearSites();
    }

    m_NX = boxSize(0);
    m_NY = boxSize(1);
    m_NZ = boxSize(2);


    m_N = boxSize;

    if (Site::nNeighborsLimit() != UNSET_UINT && check)
    {
        if (Site::nNeighborsLimit() >= min(m_N)/2)
        {
            cerr << "Neighbor reach must be lower than half the minimum box dimension to avoid sites directly affecting themselves." << endl;
            KMCSolver::exit();
        }
    }


    initializeSites();

    Site::initializeBoundaries();

    //    initializeDiffusionReactions();

    m_mainLattice->setTopology({0, m_NX,
                                0, m_NY,
                                0, m_NZ});

}

void KMCSolver::setRNGSeed(uint seedState, int defaultSeed)
{

    //    seed_type prevSeed = Seed::initialSeed;

    seed_type seed = -1;

    switch (seedState)
    {
    case Seed::specific:
        seed = static_cast<seed_type>(defaultSeed);
        break;
    case Seed::fromTime:
        seed = static_cast<seed_type>(time(NULL));
        break;
    default:
        throw std::runtime_error("Seed not specified.");
        break;
    }

    KMC_INIT_RNG(seed);

    //    if (prevSeed != Seed::initialSeed)
    //    {
    //        cout << " -- new seed set: " << seed << endl;
    //    }
}


uint KMCSolver::refCounter = 0;
