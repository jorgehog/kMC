#include "kmcsolver.h"

#include "soluteparticle.h"

#include "reactions/reaction.h"
#include "reactions/diffusion/tstdiffusion.h"
#include "reactions/diffusion/arrheniusdiffusion.h"

#include "boundary/boundary.h"

#include "ignisinterface/kmcevent.h"
#include "ignisinterface/solverevent.h"
#include "ignisinterface/kmcparticles.h"

#include <lammpswriter/lammpswriter.h>

#ifndef KMC_NO_OMP
#include <omp.h>
#endif

#include <sys/time.h>

#include <armadillo>

#include <iostream>
#include <iomanip>
#include <fstream>

#include <boost/algorithm/string.hpp>
#include <regex>

#include <numeric>

using namespace arma;
using namespace std;
using namespace kMC;

KMCSolver::KMCSolver(const Setting & root)
{

    const Setting & SolverSettings = getSetting(root, "Solver");
    const Setting & diffusionSettings = getSetting(root, {"Reactions", "Diffusion"});
    const Setting & SystemSettings = getSetting(root, "System");

    string path;

    try
    {
        path = (string)(SystemSettings["path"].c_str());
    }

    catch(const libconfig::SettingNotFoundException & exc)
    {
        path = "outfiles";
    }

    setFilepath(path);

    onConstruct();

    Reaction::loadConfig(getSetting(root, "Reactions"));

    DiffusionReaction::loadConfig(diffusionSettings);

    Site::loadConfig(SystemSettings);



    setNumberOfCycles(
                getSetting<uint>(SolverSettings, "nCycles"));

    setCyclesPerOutput(
                getSetting<uint>(SolverSettings, "cyclesPerOutput"));

    setRNGSeed(getSetting<uint>(SolverSettings, "seedType"),
               getSetting<int>(SolverSettings, "specificSeed"));

    setTargetConcentration(
                getSetting<double>(SystemSettings, "SaturationLevel"));


    uvec3 boxSize;

    boxSize(0) = getSetting(SystemSettings, "BoxSize")[0];
    boxSize(1) = getSetting(SystemSettings, "BoxSize")[1];
    boxSize(2) = getSetting(SystemSettings, "BoxSize")[2];

    setBoxSize(boxSize);

    Site::initializeBoundaries();

    initializeSites();


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


    if (m_diffusionType == DiffusionTypes::TST)
    {
        TSTDiffusion::clearAll();

    }

    else if (m_diffusionType == DiffusionTypes::Arrhenius)
    {
        ArrheniusDiffusion::clearAll();
    }


    Boundary::clearAll();


    checkAllRefCounters();


    delete m_lammpswriter;

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
    if (refCounter != 0)
    {
        exit("Instance already active.");
    }

    m_instance = this;


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

    Lattice::setCurrentParticles(new KMCParticles(this));

    m_mainLattice = new Lattice();
    m_mainLattice->setOutputPath(filePath());

    m_solverEvent = new SolverEvent();
    m_mainLattice->addEvent(m_solverEvent);

    m_lammpswriter = new lammpswriter(7, "kMC", m_filepath);

    if (m_dumpLAMMPS)
    {
        m_lammpswriter->setSystemSize(m_NX, m_NY, m_NZ);

        m_dumpFileEvent = new DumpLAMMPS();
        m_mainLattice->addEvent(m_dumpFileEvent);
    }

    else if (m_dumpXYZ)
    {
        m_dumpFileEvent = new DumpXYZ();
        m_mainLattice->addEvent(m_dumpFileEvent);
    }
}

void KMCSolver::finalizeObject()
{
    checkRefCounter();

    clearParticles();

    Site::finalizeBoundaries();

    checkAllRefCounters();


    delete m_mainLattice;

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


void KMCSolver::registerReactionChange(Reaction *reaction, const double &newRate)
{

    const double & prevRate = reaction->rate();

    if (prevRate == newRate)
    {
        BADAssBool(!((prevRate == Reaction::UNSET_RATE) && (prevRate != Reaction::UNSET_RATE)));

        return;
    }

    else if (prevRate == Reaction::UNSET_RATE)
    {

        BADAssBool(reaction->isAllowed(), "illegal reaction assigned rate.", [&reaction] ()
        {
            cout << reaction->info() << endl;
        });

        BADAss(newRate, !=, Reaction::UNSET_RATE);

        m_kTot += newRate;

        //If there is a vacancy, we simply fill it.
        if (!m_availableReactionSlots.empty())
        {

            const uint slot = m_availableReactionSlots.at(m_availableReactionSlots.size() - 1);

            reaction->setAddress(slot);

            m_allPossibleReactions.at(slot) = reaction;

            updateAccuAllRateElements(slot, m_accuAllRates.size(), newRate);

            m_availableReactionSlots.pop_back();

        }

        //if not, we make a new element
        else
        {
            m_allPossibleReactions.push_back(reaction);

            if (m_useLocalUpdating)
            {
                m_accuAllRates.push_back(m_kTot);
            }

            reaction->setAddress(m_allPossibleReactions.size()-1);

        }

    }

    else if (newRate == Reaction::UNSET_RATE)
    {

        m_kTot -= prevRate;

        BADAss(reaction->address(), !=, Reaction::UNSET_ADDRESS);
        BADAssBool(!isEmptyAddress(reaction->address()), "address is already set as empty.");

        m_availableReactionSlots.push_back(reaction->address());

        updateAccuAllRateElements(reaction->address(), m_accuAllRates.size(), -prevRate);

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

    const uint nVacancies = m_availableReactionSlots.size();

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

    BADAssBool(isEmptyAddress(dest),  "destination should be empty.");
    BADAssBool(!isEmptyAddress(orig), "origin should not be empty.");

    Reaction * swappedReaction = m_allPossibleReactions.at(orig);

    BADAss(orig,                    ==, swappedReaction->address(), "mismatch in address.", KMCBAI( swappedReaction->info()));
    BADAssBool(swappedReaction->isAllowed(), "swapped reaction should be allowed and active.", KMCBAI( swappedReaction->info()));
    BADAss(swappedReaction->rate(), !=, Reaction::UNSET_RATE, "Reaction should not appear.");

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

}

void KMCSolver::updateAccuAllRateElements(const uint from, const uint to, const double value)
{
    if (!m_useLocalUpdating)
    {
        return;
    }

    BADAss(from, <=, to);

#ifndef KMC_NO_OMP
#pragma omp parallel for
#endif
    for (uint i = from; i < to; ++i)
    {
        m_accuAllRates[i] += value;

        BADAss(m_accuAllRates.at(i), >=, -minRateThreshold(), "should be zero", KMCBAI(SoluteParticle::nParticles()));
    }

}

void KMCSolver::remakeAccuAllRates()
{
    m_kTot = 0;
    uint i = 0;

    m_accuAllRates.resize(m_allPossibleReactions.size());

    for (Reaction *r : m_allPossibleReactions)
    {
        BADAssBool(!r->hasVacantStatus(), "Reaction should be enabled.", [&] ()
        {
            cout << r->info() << endl;
        });

        BADAss(r->rate(), !=, Reaction::UNSET_RATE, "UNSET RATE IN ACTIVE REACTION.", [&] ()
        {
            cout << r->info() << endl;
        });

        BADAss(r->rate(), >, 0, "Rate should be positive.", [&] ()
        {
            cout << r->info() << endl;
        });

        m_kTot += r->rate();
        m_accuAllRates[i] = m_kTot;
        ++i;
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
    o.open(m_filepath + s.str());

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

void KMCSolver::dumpLAMMPS(const uint n)
{

    m_lammpswriter->initializeNewFile(n, SoluteParticle::nParticles());

    for (SoluteParticle *particle : m_particles)
    {
        (*m_lammpswriter) << particle->particleState()
                          << particle->x()
                          << particle->y()
                          << particle->z()
                          << particle->nNeighbors()
                          << particle->energy()
                          << particle->species();
    }

    m_lammpswriter->finalize();

}



void KMCSolver::forEachSiteDo(function<void (uint x, uint y, uint z, Site *)> applyFunction) const
{
    for (uint x = 0; x < m_NX; ++x)
    {
        for (uint y = 0; y < m_NY; ++y)
        {
            for (uint z = 0; z < m_NZ; ++z)
            {
                applyFunction(x, y, z, getSite(x, y, z));
            }
        }
    }
}


void KMCSolver::initializeSites()
{

    BADAss(Site::_refCount(), ==, 0, "Sites was not cleared properly.");

    uint xTrans, yTrans, zTrans, m_NX_full, m_NY_full, m_NZ_full;

    m_boundaryPadding = Site::nNeighborsLimit();

    m_NX_full = 2*m_boundaryPadding + m_NX;
    m_NY_full = 2*m_boundaryPadding + m_NY;
    m_NZ_full = 2*m_boundaryPadding + m_NZ;


    sites = new Site***[m_NX_full];

    for (uint x = 0; x < m_NX_full; ++x)
    {

        sites[x] = new Site**[m_NY_full];

        for (uint y = 0; y < m_NY_full; ++y)
        {

            sites[x][y] = new Site*[m_NZ_full];

            for (uint z = 0; z < m_NZ_full; ++z)
            {
                if ((x >= m_boundaryPadding && x < m_NX + m_boundaryPadding) &&
                        (y >= m_boundaryPadding && y < m_NY + m_boundaryPadding) &&
                        (z >= m_boundaryPadding && z < m_NZ + m_boundaryPadding))
                {
                    //renormalize so that Site::nNeighborsLimit() points to site 0 and so on.
                    sites[x][y][z] = new Site();
                }
            }
        }
    }


    BADAss(Site::_refCount(), !=, 0, "Can't simulate an empty system.");
    BADAss(Site::_refCount(), ==, NX()*NY()*NZ(), "Wrong number of sites initialized.");

    //Boundaries
    for (uint x = 0; x < m_NX_full; ++x)
    {
        for (uint y = 0; y < m_NY_full; ++y)
        {
            for (uint z = 0; z < m_NZ_full; ++z)
            {
                if (!((x >= m_boundaryPadding && x < m_NX + m_boundaryPadding) &&
                      (y >= m_boundaryPadding && y < m_NY + m_boundaryPadding) &&
                      (z >= m_boundaryPadding && z < m_NZ + m_boundaryPadding)))
                {


                    Boundary::setupCurrentBoundaries(x, y, z, m_boundaryPadding);

                    xTrans = Site::boundaries(0, 0)->transformCoordinate((int)x - (int)m_boundaryPadding);

                    yTrans = Site::boundaries(1, 0)->transformCoordinate((int)y - (int)m_boundaryPadding);

                    zTrans = Site::boundaries(2, 0)->transformCoordinate((int)z - (int)m_boundaryPadding);

                    if (Boundary::isBlocked(xTrans, yTrans, zTrans))
                    {
                        sites[x][y][z] = NULL;
                    }

                    else
                    {
                        sites[x][y][z] = getSite(xTrans, yTrans, zTrans);
                    }



                }
            }
        }
    }

    //Whatever fancy geometry is set up by boundaries (spheres etc.)
    uint x, y, z;
    for (uint d = 0; d < 3; ++d)
    {
        for (uint o = 0; o < 2; ++o)
        {
            const Boundary *boundary = Site::boundaries(d, o);

            if (boundary->type() != Boundary::SphericalEdge)
            {
                continue;
            }

            for (uint n = 0; n < boundary->boundarySize(); ++n)
            {
                boundary->getBoundarySite(n, x, y, z);

                boundary->applyBoundaryTransform(x, y, z, [this, boundary] (uint h, uint l, uint w)
                {
                   uint outsideBoundary = h;

                   while (outsideBoundary != boundary->bigbound())
                   {
                       outsideBoundary += boundary->orientationAsSign();

                       boundary->applyInverseBoundaryTransform(outsideBoundary, l, w, [this] (uint x, uint y, uint z)
                       {
                           sites[x + m_boundaryPadding][y + m_boundaryPadding][z + m_boundaryPadding] = NULL;
                       });
                   }
                });
            }
        }
    }

    initializeParticles();

}


void KMCSolver::clearSites()
{

    BADAss(Site::_refCount(), !=, 0, "Sites already cleared.");

    BADAss(SoluteParticle::nParticles(), ==, 0, "Cannot clear sites with particles active.");

    KMCDebugger_SetEnabledTo(false);


    uint m_NX_full = 2*m_boundaryPadding + m_NX;
    uint m_NY_full = 2*m_boundaryPadding + m_NY;
    uint m_NZ_full = 2*m_boundaryPadding + m_NZ;


    for (uint x = 0; x < m_NX_full; ++x)
    {
        for (uint y = 0; y < m_NY_full; ++y)
        {
            for (uint z = 0; z < m_NZ_full; ++z)
            {
                if ((x >= m_boundaryPadding && x < m_NX + m_boundaryPadding) &&
                        (y >= m_boundaryPadding && y < m_NY + m_boundaryPadding) &&
                        (z >= m_boundaryPadding && z < m_NZ + m_boundaryPadding))
                {
                    delete sites[x][y][z];
                }
            }

            delete [] sites[x][y];
        }

        delete [] sites[x];
    }


    delete [] sites;

    BADAss(Site::_refCount(), ==, 0, "Sites was not cleared properly.");

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

    BADAssClose(kTot, m_kTot, 1E-15);

}

uint KMCSolver::binarySearchForInterval(const double target, const vector<double> &intervals)
{

    BADAss(intervals.size(), !=, 0, "Number of intervals cannot be zero.");

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

void KMCSolver::swapMainSolverEventWith(KMCEvent *event)
{
    m_mainLattice->removeEvent(m_solverEvent->meshAddress());

    event->setManualPriority(0);
    addEvent(event);
}

double KMCSolver::getBruteForceTotalEnergy() const
{
    double eTot = 0;

    for (SoluteParticle *particle : m_particles)
    {
        eTot += particle->getBruteForceEnergy();
    }

    return eTot;
}



bool KMCSolver::spawnParticle(SoluteParticle *particle, const uint x, const uint y, const uint z, bool checkIfLegal)
{

    particle->trySite(x, y, z);

    if (particle->site() == NULL)
    {
        particle->resetSite();

        return false;
    }

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

    m_particles.push_back(particle);

    particle->setSite(x, y, z);

    BADAssBool(!checkIfLegal || particle->nNeighbors() == 0);

    return true;

}

SoluteParticle *KMCSolver::forceSpawnParticle(const uint x, const uint y, const uint z, const uint species, const bool sticky)
{
    BADAssBool(!getSite(x, y, z)->isActive());

    SoluteParticle *particle = new SoluteParticle(species, sticky);

    if (!spawnParticle(particle, x, y, z, false))
    {
        delete particle;
        return NULL;
    }

    return particle;

}

void KMCSolver::despawnParticle(SoluteParticle *particle)
{
    BADAss(particle, !=, NULL, "particle does not exist.");
    BADAss(particle->site(), !=, NULL, "despawning particle that was never initialized. (Despawning inside particle loop?)");
    BADAssBool(particle->site()->isActive(), "this should never happen.");

    SoluteParticle::popAffectedParticle(particle);

    BADAssBool(isRegisteredParticle(particle));

    m_particles.erase(std::find(m_particles.begin(), m_particles.end(), particle));

    BADAssBool(!isRegisteredParticle(particle));

    delete particle;

}

const Reaction *KMCSolver::executeRandomReaction()
{

    uint which = KMC_RNG_UNIFORM()*m_allPossibleReactions.size();

    std::sort(m_availableReactionSlots.begin(), m_availableReactionSlots.end());

    for (const uint & vacancy : m_availableReactionSlots)
    {
        if (vacancy == which)
        {
            which++;
        }
    }

    BADAss(which, <, m_allPossibleReactions.size());

    Reaction *reaction = m_allPossibleReactions.at(which);

    reaction->execute();

    return reaction;
}


void KMCSolver::initializeCrystal(const double relativeSeedSize, const uint species, const bool sticky)
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
                            forceSpawnParticle(i, j, k, species, sticky);
                        }
                    }
                }

            }
        }
    }

    KMCDebugger_ResetEnabled();

}

void KMCSolver::initializeSolutionBath(const uint species, const bool sticky)
{
    const double margin = 1.75;
    uint effectiveVolume = 8; //eV = (difflength + 1)^3

    effectiveVolume *= margin;

    uint NFree = NX()*NY()*NZ() - SoluteParticle::nParticles();

    uint N = NFree*targetConcentration();

    if (N > NFree/effectiveVolume)
    {
        cerr << "Not enough space to place " << N << "particles sufficiently apart from eachother with concentration " << m_targetConcentration << " Maximum concentration: " << 1./effectiveVolume << endl;
        exit();
    }


    insertRandomParticles(N, species, sticky, true);

}

void KMCSolver::initializeLayers(const uint height, const uint start, const uint species, const bool sticky)
{
    for (uint x = 0; x < m_NX; ++x)
    {
        for (uint y = 0; y < m_NY; ++y)
        {
            for (uint h = 0; h < height; ++h)
            {
                forceSpawnParticle(x, y, start + h, species, sticky);
            }
        }
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
        BADAssEqual(t.at(i), p->particleStateShortName());
        BADAssEqual(p->nNeighborsSum(), nn.at(i));
        BADAssClose(p->energy(), e.at(i), 1E-5);
        ++i;
    }
#endif

    if (m_dumpXYZ || m_dumpLAMMPS)
    {
        m_dumpFileEvent->setOffset(frame + 1);
    }

}

void KMCSolver::initializeFromLAMMPS(string path, uint frame)
{
    const string &prevPath = m_lammpswriter->path();

    m_lammpswriter->setPath(path);
    m_lammpswriter->loadFile(frame);

    clearSites();
    setBoxSize(m_lammpswriter->systemSizeX(),
               m_lammpswriter->systemSizeY(),
               m_lammpswriter->systemSizeZ());
    initializeSites();


    double ePot = 0;
    uint type = 0, species = 0, nNeighbors = 0, x = 0, y = 0, z = 0;

    for (uint i = 0; i < m_lammpswriter->nParticles(); ++i)
    {
        (*m_lammpswriter) >> type
                          >> x
                          >> y
                          >> z
                          >> nNeighbors
                          >> ePot
                          >> species;

        forceSpawnParticle(x, y, z, species);
    }

    m_lammpswriter->finalize();

    m_lammpswriter->setPath(prevPath);

    if (m_dumpXYZ || m_dumpLAMMPS)
    {
        m_dumpFileEvent->setOffset(frame + 1);
    }

}

void KMCSolver::insertRandomParticles(const uint N, const uint species, const bool sticky, bool checkIfLegal)
{
    uint n = 0;
    while (n != N)
    {
        insertRandomParticle(species, sticky, checkIfLegal);
        n++;
    }
}

void KMCSolver::insertRandomParticle(const uint species, const bool sticky, bool checkIfLegal)
{

    const double margin = 1.75;
    uint effectiveVolume = 8; //eV = (difflength + 1)^3

    effectiveVolume *= margin;

    uint NFree = NX()*NY()*NZ() - SoluteParticle::nParticles();

    uint N = NFree*targetConcentration();

    if (N > NFree/effectiveVolume)
    {
        cerr << "Not enough space to place " << N << "particles sufficiently apart from eachother with concentration " << m_targetConcentration << " Maximum concentration: " << 1./effectiveVolume << endl;
        exit();
    }

    SoluteParticle *particle = new SoluteParticle(species, sticky);

    forceRandomPosition(particle, checkIfLegal);

}

void KMCSolver::forceRandomPosition(SoluteParticle *particle, bool checkIfLegal)
{

    uint N = 10000;

    uint x, y, z;
    bool spawned = false;

    uint c = 0;
    while (!spawned)
    {
        x = KMC_RNG_UNIFORM()*NX();

        y = KMC_RNG_UNIFORM()*NY();

        z = KMC_RNG_UNIFORM()*NZ();

        spawned = spawnParticle(particle, x, y, z, checkIfLegal);

        c++;

        if (c == N)
        {
            exit("Unable to place particle.");
        }

    }

}

void KMCSolver::rotateSystem(const double yaw, const double pitch, const double roll)
{

    using std::sin;
    using std::cos;

    //Convert to radians
    const double y = yaw/180.0*datum::pi;
    const double p = pitch/180.0*datum::pi;
    const double r = roll/180.0*datum::pi;

    //Geometric values
    const double sy = sin(y);
    const double sp = sin(p);
    const double sr = sin(r);

    const double cy = cos(y);
    const double cp = cos(p);
    const double cr = cos(r);


    //Rotational matrices
    mat Y(3, 3, fill::eye);
    mat P(3, 3, fill::eye);
    mat R(3, 3, fill::eye);

    Y(1, 1) =  cy; Y(1, 2) = sy;
    Y(2, 1) = -sy; Y(2, 2) = cy;

    P(0, 0) =  cp; P(0, 2) = sp;
    P(2, 0) = -sp; P(2, 2) = cp;

    R(0, 0) =  cr; R(0, 1) = sr;
    R(1, 0) = -sr; R(1, 1) = cr;

    cout << Y << P << R << endl;


    //Total transformation matrix
    mat RotMat = Y*P*R;

    cout << RotMat << endl;

    mat transformedPositions(3, SoluteParticle::nParticles());

    uint c = 0;
    for (SoluteParticle *particle : m_particles)
    {
        ivec3 centerReferencePos = {(int)particle->x() - (int)NX()/2,
                                    (int)particle->y() - (int)NY()/2,
                                    (int)particle->z() - (int)NZ()/2};

        cout << centerReferencePos.t() << endl;

        vec3 transPos = RotMat*centerReferencePos;

        cout << transPos.t() << "\n--------------" << endl;

        particle->disableSite();

        transformedPositions.col(c) = transPos;

        c++;
    }

    cout << "---" << endl;

    vector<ivec3> collisions;
    vector<SoluteParticle*> colliders;

    c = 0;
    for (SoluteParticle *particle : m_particles)
    {
        int tx = std::round(transformedPositions(0, c) + NX()/2);
        int ty = std::round(transformedPositions(1, c) + NY()/2);
        int tz = std::round(transformedPositions(2, c) + NZ()/2);

        cout << transformedPositions.col(c).t();
        cout << tx << " " << ty << " " << tz;

        if (getSite(tx, ty, tz)->isActive())
        {
            cout << " " << "collided.";
            collisions.push_back({tx, ty, tz});
            colliders.push_back(particle);
        }

        else
        {
            particle->setSite(tx, ty, tz);
        }

        cout << endl;

        dumpLAMMPS(1337 + c + 1);

        c++;
    }

    int xd, yd, zd;
    for (uint i = 0; i < colliders.size(); ++i)
    {
        int tx = collisions.at(i)(0);
        int ty = collisions.at(i)(1);
        int tz = collisions.at(i)(2);


        double min = numeric_limits<double>::max();

        Site::forShellDo(tx, ty, tz, 1, [&] (Site *site, int dx, int dy, int dz)
        {
            if (site->isActive())
            {
                return;
            }

            int x = Boundary::currentBoundaries(0)->transformCoordinate(tx + dx);
            int y = Boundary::currentBoundaries(1)->transformCoordinate(ty + dy);
            int z = Boundary::currentBoundaries(2)->transformCoordinate(tz + dz);

            double dr2 = (tx - x)*(tx - x) + (ty - y)*(ty - y) + (tz - z)*(tz - z);

            if (dr2 < min)
            {
                min = dr2;

                xd = x;
                yd = y;
                zd = z;

            }

        });

        colliders.at(i)->setSite(0, 0, i);

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

    if (!m_useLocalUpdating)
    {
        remakeAccuAllRates();
    }

}

const uint &KMCSolver::cycle() const
{
    return m_solverEvent->cycle();
}


void KMCSolver::setBoxSize(const uint NX, const uint NY, const uint NZ, bool check)
{

    BADAss(Site::_refCount(), ==, 0, "Sites need to be c.leared before a new boxsize is set.");

    m_NX = NX;
    m_NY = NY;
    m_NZ = NZ;

    m_N = {NX, NY, NZ};

    if (Site::nNeighborsLimit() != UNSET_UINT && check)
    {
        for (uint i = 0; i < 3; ++i)
        {
            if (Site::boundaries(i, 0)->type() == Boundary::BoundaryTypes::Periodic)
            {
                if (Site::nNeighborsLimit() >= m_N(i)/2)
                {
                    cerr << "Neighbor reach must be lower than half the box dimension (periodic only) to avoid sites directly affecting themselves." << endl;
                    KMCSolver::exit();
                }
            }
        }
    }

    m_mainLattice->setTopology({0, 0, 0, NX, NY, NZ});

    m_lammpswriter->setSystemSize(NX, NY, NZ);

}

void KMCSolver::resetBoxSize(const uint NX, const uint NY, const uint NZ, bool check)
{
    Site::finalizeBoundaries();

    clearSites();
    setBoxSize(NX, NY, NZ, check);
    initializeSites();

    Site::initializeBoundaries();
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
        particle = NULL;
    }

    if (m_useLocalUpdating)
    {
        BADAssClose(0, std::accumulate(m_accuAllRates.begin(), m_accuAllRates.end(), 0.0), 1E-5);
    }

    SoluteParticle::clearAll();

    m_particles.clear();

    BADAss(accu(SoluteParticle::totalParticlesVector()), ==, 0);
    BADAss(SoluteParticle::nParticles(), ==, 0);
    BADAss(SoluteParticle::affectedParticles().size(), ==, 0);
    BADAssClose(SoluteParticle::totalEnergy(), 0, 1E-5);

    SoluteParticle::setZeroTotalEnergy();

    m_allPossibleReactions.clear();

    m_accuAllRates.clear();

    m_availableReactionSlots.clear();

}

KMCSolver *KMCSolver::m_instance = NULL;

bool KMCSolver::m_useLocalUpdating = true;
bool KMCSolver::m_dumpXYZ = false;
bool KMCSolver::m_dumpLAMMPS = true;

uint KMCSolver::refCounter = 0;


