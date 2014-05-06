#include "kmcsolver.h"

#include "soluteparticle.h"

#include "reactions/reaction.h"
#include "reactions/diffusion/diffusionreaction.h"

#include "boundary/boundary.h"

#include "ignisinterface/kmcevent.h"
#include "ignisinterface/solverevent.h"
#include "ignisinterface/kmcparticles.h"

#include <omp.h>

#include <sys/time.h>

#include <armadillo>

#include <iostream>
#include <iomanip>
#include <fstream>

#include <boost/algorithm/string.hpp>
#include <regex>

using namespace arma;
using namespace std;
using namespace kMC;

KMCSolver::KMCSolver(const Setting & root)
{

    onConstruct();

    const Setting & SystemSettings = getSetting(root, "System");
    const Setting & SolverSettings = getSetting(root, "Solver");
    const Setting & diffusionSettings = getSetting(root, {"Reactions", "Diffusion"});


    Reaction::loadConfig(getSetting(root, "Reactions"));

    DiffusionReaction::loadConfig(diffusionSettings);

    Site::loadConfig(SystemSettings);


    setNumberOfCycles(
                getSetting<uint>(SolverSettings, "nCycles"));

    setCyclesPerOutput(
                getSetting<uint>(SolverSettings, "cyclesPerOutput"));

    setRNGSeed(
                getSetting<uint>(SolverSettings, "seedType"),
                getSetting<int>(SolverSettings, "specificSeed"));

    setTargetConcentration(
                getSetting<double>(SystemSettings, "SaturationLevel"));


    uvec3 boxSize;

    boxSize(0) = getSetting(SystemSettings, "BoxSize")[0];
    boxSize(1) = getSetting(SystemSettings, "BoxSize")[1];
    boxSize(2) = getSetting(SystemSettings, "BoxSize")[2];

    setBoxSize(boxSize);

    initializeSites();

    Site::initializeBoundaries();

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


    checkAllRefCounters();

    refCounter--;

}

void KMCSolver::reset()
{

    finalizeObject();


    KMCDebugger_Init();

    setRNGSeed(Seed::specific, Seed::initialSeed);

    Site::initializeBoundaries();

    setupMainLattice();

    m_kTot = 0;

}

void KMCSolver::onConstruct()
{

    m_NX = UNSET_UINT;
    m_NY = UNSET_UINT;
    m_NZ = UNSET_UINT;

    m_targetConcentration = 0;

    m_kTot = 0;

    Boundary::setMainSolver(this);

    Reaction::setMainSolver(this);

    Site::setMainSolver(this);

    SoluteParticle::setMainSolver(this);


    setupMainLattice();


    refCounter++;

}

void KMCSolver::setupMainLattice()
{

    MainLattice::setCurrentParticles(new KMCParticles(this));

    m_mainLattice = new MainLattice();

    m_solverEvent = new SolverEvent();
    m_mainLattice->addEvent(m_solverEvent);

    if (m_dumpXYZ)
    {
        xyzEvent = new DumpXYZ();
        m_mainLattice->addEvent(xyzEvent);
    }
}

void KMCSolver::finalizeObject()
{
    checkRefCounter();

    clearParticles();


    Site::finalizeBoundaries();

    checkAllRefCounters();


    delete m_mainLattice;

    Event<uint>::resetEventParameters();


    KMCDebugger_Finalize();
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

void KMCSolver::checkAllRefCounters()
{
    if (SoluteParticle::nParticles() != 0)
    {
        stringstream s;
        s << "After deletion: " << SoluteParticle::nParticles() << " particles active.";
        exit(s.str());
    }


    if (Reaction::_refCount() != 0)
    {
        stringstream s;
        s << "After deletion: " << Reaction::_refCount() << " reactions active.";
        exit(s.str());
    }

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

    KMCDebugger_Assert(orig,                    ==, swappedReaction->address(), "mismatch in address.", swappedReaction->getFinalizingDebugMessage());
    KMCDebugger_AssertBool(swappedReaction->isAllowed(), "swapped reaction should be allowed and active.", swappedReaction->getFinalizingDebugMessage());

    m_allPossibleReactions.at(dest) = swappedReaction;

    swappedReaction->setAddress(dest);

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

void KMCSolver::updateAccuAllRateElements(const uint from, const uint to, const double value)
{
    KMCDebugger_Assert(from, <=, to);

#ifndef KMC_NO_OMP
#pragma omp parallel for
#endif
    for (uint i = from; i < to; ++i)
    {
        m_accuAllRates[i] += value;

        KMCDebugger_Assert(m_accuAllRates.at(i), >=, -minRateThreshold());
    }

}


bool KMCSolver::isEmptyAddress(const uint address) const
{
    return std::find(m_availableReactionSlots.begin(), m_availableReactionSlots.end(), address)
            != m_availableReactionSlots.end();
}

bool KMCSolver::isRegisteredParticle(SoluteParticle *particle) const
{
    return std::find(m_particles.begin(), m_particles.end(), particle) != m_particles.end();
}

bool KMCSolver::isPossibleReaction(Reaction *reaction) const
{
    return std::find(m_allPossibleReactions.begin(), m_allPossibleReactions.end(), reaction) != m_allPossibleReactions.end();
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

void KMCSolver::dumpXYZ(const uint n)
{

    stringstream s;
    s << "kMC" << n << ".xyz";

    ofstream o;
    o.open("outfiles/" + s.str());



    stringstream surface;
    stringstream crystal;
    stringstream solution;

    s.str(string());

    for (SoluteParticle *particle : particles())
    {

        s << "\n"
          << particle->particleStateShortName() << " "
          << particle->x() << " " << particle->y() << " " << particle->z() << " "
          << particle->nNeighbors() << " "
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

    }

    o << particles().size() << "\n";
    o << m_NX << " " << m_NY << " " << m_NZ;
    o << surface.str() << crystal.str() << solution.str();
    o.close();

}


void KMCSolver::forEachSiteDo(function<void (uint x, uint y, uint z, Site *)> applyFunction) const
{
    for (uint x = 0; x < m_NX; ++x)
    {
        for (uint y = 0; y < m_NY; ++y)
        {
            for (uint z = 0; z < m_NZ; ++z)
            {
                KMCDebugger_Assert(getSite(x, y, z), !=, NULL);
                applyFunction(x, y, z, getSite(x, y, z));
            }
        }
    }
}


void KMCSolver::initializeSites()
{

    KMCDebugger_Assert(Site::_refCount(), ==, 0, "Sites was not cleared properly.");

    uint xTrans, yTrans, zTrans, m_NX_full, m_NY_full, m_NZ_full;


    m_NX_full = 2*Site::nNeighborsLimit() + m_NX;
    m_NY_full = 2*Site::nNeighborsLimit() + m_NY;
    m_NZ_full = 2*Site::nNeighborsLimit() + m_NZ;


    sites = new Site***[m_NX_full];

    for (uint x = 0; x < m_NX_full; ++x)
    {

        sites[x] = new Site**[m_NY_full];

        for (uint y = 0; y < m_NY_full; ++y)
        {

            sites[x][y] = new Site*[m_NZ_full];

            for (uint z = 0; z < m_NZ_full; ++z)
            {
                if ((x >= Site::nNeighborsLimit() && x < m_NX + Site::nNeighborsLimit()) &&
                        (y >= Site::nNeighborsLimit() && y < m_NY + Site::nNeighborsLimit()) &&
                        (z >= Site::nNeighborsLimit() && z < m_NZ + Site::nNeighborsLimit()))
                {
                    //renormalize so that Site::nNeighborsLimit() points to site 0 and so on.
                    sites[x][y][z] = new Site();
                }
            }
        }
    }


    KMCDebugger_Assert(Site::_refCount(), !=, 0, "Can't simulate an empty system.");
    KMCDebugger_Assert(Site::_refCount(), ==, NX()*NY()*NZ(), "Wrong number of sites initialized.");

    //Boundaries
    for (uint x = 0; x < m_NX_full; ++x)
    {
        for (uint y = 0; y < m_NY_full; ++y)
        {
            for (uint z = 0; z < m_NZ_full; ++z)
            {
                if (!((x >= Site::nNeighborsLimit() && x < m_NX + Site::nNeighborsLimit()) &&
                      (y >= Site::nNeighborsLimit() && y < m_NY + Site::nNeighborsLimit()) &&
                      (z >= Site::nNeighborsLimit() && z < m_NZ + Site::nNeighborsLimit())))
                {


                    Boundary::setupCurrentBoundaries(x, y, z, Site::nNeighborsLimit());

                    xTrans = Site::boundaries(0, 0)->transformCoordinate((int)x - (int)Site::nNeighborsLimit());

                    yTrans = Site::boundaries(1, 0)->transformCoordinate((int)y - (int)Site::nNeighborsLimit());

                    zTrans = Site::boundaries(2, 0)->transformCoordinate((int)z - (int)Site::nNeighborsLimit());

                    if (Boundary::isBlocked(xTrans, yTrans, zTrans))
                    {
                        sites[x][y][z] = NULL;
                    }

                    else
                    {
                        sites[x][y][z] = sites[xTrans + Site::nNeighborsLimit()]
                                [yTrans + Site::nNeighborsLimit()]
                                [zTrans + Site::nNeighborsLimit()];
                    }



                }
            }
        }
    }

    initializeParticles();

}


void KMCSolver::clearSites()
{

    KMCDebugger_Assert(Site::_refCount(), !=, 0, "Sites already cleared.");

    KMCDebugger_Assert(SoluteParticle::nParticles(), ==, 0, "Cannot clear sites with particles active.");

    KMCDebugger_SetEnabledTo(false);


    uint m_NX_full = 2*Site::nNeighborsLimit() + m_NX;
    uint m_NY_full = 2*Site::nNeighborsLimit() + m_NY;
    uint m_NZ_full = 2*Site::nNeighborsLimit() + m_NZ;


    for (uint x = 0; x < m_NX_full; ++x)
    {
        for (uint y = 0; y < m_NY_full; ++y)
        {
            for (uint z = 0; z < m_NZ_full; ++z)
            {
                if ((x >= Site::nNeighborsLimit() && x < m_NX + Site::nNeighborsLimit()) &&
                        (y >= Site::nNeighborsLimit() && y < m_NY + Site::nNeighborsLimit()) &&
                        (z >= Site::nNeighborsLimit() && z < m_NZ + Site::nNeighborsLimit()))
                {
                    delete sites[x][y][z];
                }
            }

            delete [] sites[x][y];
        }

        delete [] sites[x];
    }


    delete [] sites;

    KMCDebugger_Assert(Site::_refCount(), ==, 0, "Sites was not cleared properly.");

    KMCDebugger_ResetEnabled();

}

void KMCSolver::resetLastReaction()
{
    m_solverEvent->resetReaction();
}

void KMCSolver::sortReactionsByRate()
{
    std::sort(m_allPossibleReactions.begin(),
              m_allPossibleReactions.end(),
              [] (const Reaction * r1, const Reaction * r2)
              {
                    return r1->rate() < r2->rate();
              });


    double kTot = 0;

    uint address = 0;

    for (Reaction *r : m_allPossibleReactions)
    {
        r->setAddress(address);

        kTot += r->rate();

        m_accuAllRates.at(address) = kTot;

        address++;
    }

    KMCDebugger_AssertClose(kTot, m_kTot, 1E-15);

}

uint KMCSolver::binarySearchForInterval(const double target, const vector<double> &intervals)
{

    KMCDebugger_Assert(intervals.size(), !=, 0, "Number of intervals cannot be zero.");

    uint imax = intervals.size() - 1;
    uint MAX = imax;

    uint imin = 0;
    uint imid;

    // continue searching while imax != imin + 1
    do
    {

        // calculate the midpoint for roughly equal partition
        imid = imin + (imax - imin)/2;

        //Is the upper limit above mid?
        if (target > intervals[imid])
        {

            //This means that the target is the last interval.
            if (imid == MAX)
            {
                return MAX;
            }

            //Are we just infront of the limit?
            else if (target < intervals[imid + 1])
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

            //This means that the target is the first inteval.
            if (imid == 0)
            {
                return 0;
            }

            imax = imid;
        }


    } while (imid != imin);

    //If we get here, imid = imin, which means that imax = imid + 1 (deduced by integer division).
    //We choose the max value as out match.
    return imid + 1;

}



bool KMCSolver::spawnParticle(SoluteParticle *particle, const uint x, const uint y, const uint z, bool checkIfLegal)
{

    particle->trySite(x, y, z);

    if (particle->site()->isActive())
    {
        particle->resetSite();

        return false;
    }

    if (checkIfLegal)
    {

        if (!particle->isLegalToSpawn())
        {

            particle->resetSite();

            return false;
        }

    }

    particle->setSite(x, y, z);

    KMCDebugger_AssertBool(!checkIfLegal || particle->nNeighbors() == 0);

    m_particles.push_back(particle);

    return true;

}

void KMCSolver::forceSpawnParticle(const uint x, const uint y, const uint z)
{
    KMCDebugger_AssertBool(!getSite(x, y, z)->isActive());

    SoluteParticle *particle = new SoluteParticle();

    spawnParticle(particle, x, y, z, false);

}

void KMCSolver::despawnParticle(SoluteParticle *particle)
{
    KMCDebugger_Assert(particle, !=, NULL, "particle does not exist.");
    KMCDebugger_AssertBool(particle->site()->isActive(), "this should never happen.");

    SoluteParticle::popAffectedParticle(particle);

    KMCDebugger_AssertBool(isRegisteredParticle(particle));

    m_particles.erase(std::find(m_particles.begin(), m_particles.end(), particle));

    KMCDebugger_AssertBool(!isRegisteredParticle(particle));

    delete particle;

}


void KMCSolver::initializeCrystal(const double relativeSeedSize)
{

    if (relativeSeedSize >= 1.0)
    {
        KMCSolver::exit("The seed size cannot be or exceed the box size.");
    }

    else if (relativeSeedSize < 0)
    {
        KMCSolver::exit("The seed size cannot be negative.");
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
                            forceSpawnParticle(i, j, k);
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

    uint x, y, z;
    bool spawned;


    const double margin = 0.5;
    uint effectiveVolume = 8; //eV = (difflength + 1)^3

    effectiveVolume *= margin;


    uint NFree = NX()*NY()*NZ() - SoluteParticle::nParticles();

    uint N = NFree*targetConcentration();

    if (N > NFree/effectiveVolume)
    {
        cerr << "Not enough space to place " << N << "particles sufficiently apart from eachother with concentration " << m_targetConcentration << " Maximum concentration: " << 1./effectiveVolume << endl;
        exit();
    }

    uint n = 0;

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

void KMCSolver::initializeFromXYZ(string path, uint frame)
{

    cout << "initializing from XYZ is highly experimental." << endl;

    string line;
    vector<string> tokens;

    ifstream file;

    stringstream filename;
    filename << "kMC" << frame << ".xyz";

    string fullpath = path + "/" + filename.str();

    file.open(fullpath);

    if (!file.good())
    {
        exit("file doest exist: " + fullpath);
    }

    getline(file, line);

    uint N = atoi(line.c_str());

    getline(file, line);

    boost::split(tokens, line, boost::is_any_of(" "));

    uvec3 boxSize;

    bool initBox = true;
    for (uint i = 0; i < 3; ++i)
    {
        boxSize(i) = atoi(tokens.at(i).c_str());

        if (boxSize(i) == 0)
        {
            cout << "Warning: Unable to deduce system size from file.";

            initBox = false;
            break;
        }
    }

    if (initBox)
    {
        cout << boxSize << endl;

        Site::finalizeBoundaries();

        clearSites();
        setBoxSize(boxSize);
        initializeSites();

        Site::initializeBoundaries();

    }

    uint x, y, z;

    vector<uint> nn;
    vector<double> e;
    vector<string> t;

    for (uint i = 0; i < N; ++i)
    {
        getline(file, line);

        boost::split(tokens, line, boost::is_any_of(" "));

        t.push_back(tokens.at(0));

        x = atoi(tokens.at(1).c_str());
        y = atoi(tokens.at(2).c_str());
        z = atoi(tokens.at(3).c_str());

        nn.push_back(atoi(tokens.at(4).c_str()));
        e.push_back(atof(tokens.at(5).c_str()));

        forceSpawnParticle(x, y, z);

    }


#ifndef KMC_NO_DEBUG
    uint i = 0;
    for (SoluteParticle *p : m_particles)
    {
        KMCDebugger_AssertEqual(t.at(i), p->particleStateShortName());
        KMCDebugger_AssertEqual(p->nNeighborsSum(), nn.at(i));
        KMCDebugger_AssertClose(p->energy(), e.at(i), DiffusionReaction::potentialBox().min());
        ++i;
    }
#endif

    if (m_dumpXYZ)
    {
        xyzEvent->setOffset(frame + 1);
    }

}


void KMCSolver::initializeParticles()
{
    for (SoluteParticle *particle : m_particles)
    {
        particle->setVectorSizes();
        particle->setupAllNeighbors();
    }
}




void KMCSolver::getRateVariables()
{

    SoluteParticle::updateAffectedParticles();

    reshuffleReactions();

    KMCDebugger_AssertClose(accuAllRates().at(0), allPossibleReactions().at(0)->rate(), minRateThreshold(), "zeroth accuallrate should be the first rate.");

}


void KMCSolver::setBoxSize(const uint NX, const uint NY, const uint NZ, bool check)
{

    KMCDebugger_Assert(Site::_refCount(), ==, 0, "Sites need to be c.leared before a new boxsize is set.");

    m_NX = NX;
    m_NY = NY;
    m_NZ = NZ;

    m_N = {NX, NY, NZ};

    if (Site::nNeighborsLimit() != UNSET_UINT && check)
    {
        if (Site::nNeighborsLimit() >= min(m_N)/2)
        {
            cerr << "Neighbor reach must be lower than half the minimum box dimension to avoid sites directly affecting themselves." << endl;
            KMCSolver::exit();
        }
    }

    m_mainLattice->setTopology({0, m_NX,
                                0, m_NY,
                                0, m_NZ});

}

void KMCSolver::setRNGSeed(uint seedState, int defaultSeed)
{

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

}

void KMCSolver::clearParticles()
{

    for (SoluteParticle *particle : m_particles)
    {
        delete particle;
    }

    SoluteParticle::clearAll();

    m_particles.clear();

    KMCDebugger_Assert(accu(SoluteParticle::totalParticlesVector()), ==, 0);
    KMCDebugger_Assert(SoluteParticle::nParticles(), ==, 0);
    KMCDebugger_Assert(SoluteParticle::affectedParticles().size(), ==, 0);
    KMCDebugger_AssertClose(SoluteParticle::totalEnergy(), 0, 1E-5);

    SoluteParticle::setZeroTotalEnergy();

    m_allPossibleReactions.clear();

    m_accuAllRates.clear();

    m_availableReactionSlots.clear();

}


bool KMCSolver::m_dumpXYZ = true;

uint KMCSolver::refCounter = 0;
