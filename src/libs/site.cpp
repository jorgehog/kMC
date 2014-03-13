#include "site.h"
#include "kmcsolver.h"
#include "reactions/reaction.h"
#include "reactions/diffusion/diffusionreaction.h"

#include "boundary/periodic/periodic.h"
#include "boundary/edge/edge.h"
#include "boundary/wall/surface.h"
#include "boundary/concentrationwall/concentrationwall.h"

#include "debugger/debugger.h"

using namespace kMC;


Site::Site(uint _x, uint _y, uint _z) :
    m_nNeighborsSum(0),
    m_active(false),
    m_isFixedCrystalSeed(false),
    m_x(_x),
    m_y(_y),
    m_z(_z),
    m_r({_x, _y, _z}),
    m_energy(0),
    m_particleState(ParticleStates::solution)
{

}


Site::~Site()
{

    clearNeighborhood();

    clearAllReactions();

    if (isActive())
    {
        m_totalActiveSites--;
    }

}




void Site::updateAffectedSites()
{

    for (Site* site : m_affectedSites)
    {
        KMCDebugger_AssertBool(site->isActive(), "Affected site should be active.", site->info());
        site->updateReactions();
        site->calculateRates();
    }

    m_affectedSites.clear();

}

void Site::selectUpdateFlags()
{
    for (Site * site : m_affectedSites)
    {

        for (Reaction * reaction : site->siteReactions())
        {
            reaction->selectTriumphingUpdateFlag();
        }

    }
}


void Site::setParticleState(int newState)
{

    KMCDebugger_MarkPre(particleStateName());


    KMCDebugger_Assert(newState, !=, m_particleState, "switching particle states to same state...", info());

    /* ################## crystal -> surface | solution ################### */
    if (m_particleState == ParticleStates::crystal)
    {

        decrystallize();

        queueAffectedSites();

    }
    /* #################################################################### */

    else if (m_particleState == ParticleStates::surface)
    {


        /* ##################     surface -> crystal      ################### */
        if (newState == ParticleStates::crystal)
        {

            KMCDebugger_AssertBool(isActive());

            crystallize();

        }
        /* ################################################################## */




        /* ##################     surface -> solution     ################### */
        else if (newState == ParticleStates::solution)
        {

            KMCDebugger_AssertBool(!isActive(), "surface should not be active.", info());

            if (!qualifiesAsSurface())
            {
                m_particleState = ParticleStates::solution;
                queueAffectedSites();
            }


            KMCDebugger_PushImplication(this, particleStateName().c_str());
        }
        /* ################################################################## */



    }

    else if (m_particleState == ParticleStates::solution)
    {


        KMCDebugger_Assert(newState, !=, ParticleStates::crystal, "should never happen.", info());
        KMCDebugger_Assert(newState, ==, ParticleStates::surface, "should allways happen.", info());


        /* ##################     solution -> surface     ################### */

        //if a particle is present, we crystallize it
        if (m_active)
        {
            crystallize();
        }

        else
        {
            m_particleState = ParticleStates::surface;

            KMCDebugger_PushImplication(this, particleStateName().c_str());

            queueAffectedSites();
        }
        /* ################################################################## */


    }

    KMCDebugger_AssertBool(!(isSurface() && isActive()), "surface should not be active.", info());

}


//All reactions must be legal if site is allowed to spawn.
bool Site::isLegalToSpawn()
{

    if (m_active)
    {
        return false;
    }

    for (Reaction * r : m_siteReactions)
    {
        if (!r->isAllowed())
        {
            return false;
        }
    }


    return true;

}

bool Site::qualifiesAsCrystal()
{

    if (isFixedCrystalSeed())
    {
        return true;
    }

    else if (hasNeighboring(ParticleStates::fixedCrystal, 1))
    {
        return true;
    }

    else if (countNeighboring(ParticleStates::crystal, 1) >= m_nNeighborsToCrystallize)
    {
        return true;
    }

    return false;

}

bool Site::qualifiesAsSurface()
{
    return !isActive() && hasNeighboring(ParticleStates::crystal, DiffusionReaction::separation());
}


void Site::loadConfig(const Setting &setting)
{

    setNNeighborsToCrystallize(getSurfaceSetting<uint>(setting, "nNeighboursToCrystallize"));

    setNNeighborsLimit(getSurfaceSetting<uint>(setting, "nNeighborsLimit"));

}

void Site::initializeBoundaries()
{
    bool enabled = KMCDebugger_IsEnabled;
    KMCDebugger_SetEnabledTo(false);

    for (uint i = 0; i < 3; ++i) {
        for (uint j = 0; j < 2; ++j) {
            m_boundaries(i, j)->initialize();
        }
    }

    KMCDebugger_SetEnabledTo(enabled);
}

void Site::updateBoundaries()
{
    for (uint i = 0; i < 3; ++i) {
        for (uint j = 0; j < 2; ++j) {
            m_boundaries(i, j)->update();
        }
    }
}


void Site::addReaction(Reaction *reaction)
{
    m_siteReactions.push_back(reaction);
}

void Site::updateReactions()
{

    m_activeReactions.clear();

    if (!m_active)
    {
        return;
    }

    for (Reaction* reaction : m_siteReactions)
    {
        if (reaction->isAllowed())
        {
            m_activeReactions.push_back(reaction);
        }
    }

}

void Site::clearAllReactions()
{
    for (Reaction * reaction : m_siteReactions)
    {
        delete reaction;
    }

    m_siteReactions.clear();

    m_activeReactions.clear();

}

void Site::spawnAsFixedCrystal()
{
    m_particleState = ParticleStates::surface;
    m_isFixedCrystalSeed = true;

    clearAllReactions();

    activate();

    m_particleState = ParticleStates::fixedCrystal;

}

void Site::deactivateFixedCrystal()
{

    KMCDebugger_Assert(particleState(), ==, ParticleStates::fixedCrystal);

    m_isFixedCrystalSeed = false;

    m_particleState = ParticleStates::crystal;

    deactivate();

    initializeDiffusionReactions();

}

void Site::crystallize()
{

    KMCDebugger_AssertBool(isActive(), "should not attempt to crystallize an unactive site.", info());

    if (qualifiesAsCrystal())
    {
        //No need to test if it has neigh crystals because
        //diffusion reactions deactivates old spot before activating new spot.
        //Which will remove access surface.
        m_particleState  = ParticleStates::crystal;

        KMCDebugger_PushImplication(this, particleStateName().c_str());
        propagateToNeighbors(ParticleStates::solution, ParticleStates::surface, DiffusionReaction::separation());
    }
    else
    {
        KMCDebugger_AssertBool(!(DiffusionReaction::separation() == 1 && Site::nNeighborsToCrystallize() == 1), "With a single neighboring crystal needed, everything should crystallize if asked.", info());

        m_particleState = ParticleStates::solution;

        KMCDebugger_PushImplication(this, particleStateName().c_str());
    }
}

void Site::decrystallize()
{

    if (qualifiesAsSurface())
    {
        m_particleState = ParticleStates::surface;
        propagateToNeighbors(ParticleStates::any, ParticleStates::solution, DiffusionReaction::separation());
    }

    else if (!qualifiesAsCrystal())
    {
        m_particleState = ParticleStates::solution;
        propagateToNeighbors(ParticleStates::any, ParticleStates::solution, DiffusionReaction::separation());
    }

    KMCDebugger_PushImplication(this, particleStateName().c_str());

}


void Site::calculateRates()
{
    for (Reaction* reaction : m_activeReactions)
    {
        reaction->selectTriumphingUpdateFlag();
        reaction->calcRate();
    }
}

void Site::initializeDiffusionReactions()
{


    KMCDebugger_Assert(m_siteReactions.size(), ==, 0, "Sitereactions are already set", info());
    KMCDebugger_AssertBool(!isActive() || isFixedCrystalSeed(), "Non FixedCrystal Site should not be active when reactions are initialized.", info());

    if (isFixedCrystalSeed())
    {
        return;
    }

    Site * destination;

    //For each site, loop over all closest neighbors
    for (uint i = 0; i < 3; ++i)
    {
        for (uint j = 0; j < 3; ++j)
        {
            for (uint k = 0; k < 3; ++k)
            {

                destination = neighborHood(Site::nNeighborsLimit() - 1 + i,
                                           Site::nNeighborsLimit() - 1 + j,
                                           Site::nNeighborsLimit() - 1 + k);

                //This means that the destination is blocked by boundaries
                if (destination != NULL)
                {
                    //This menas we are not at the current site.
                    if(destination != this)
                    {
                        DiffusionReaction* diffusionReaction = new DiffusionReaction(this, destination);
                        addReaction(diffusionReaction);
                    }

                    else
                    {
                        KMCDebugger_AssertBool((i == 1) && (j == 1) && (k == 1));
                    }
                }

            }
        }
    }
}


void Site::setMainSolver(KMCSolver *solver)
{
    m_solver = solver;
}


void Site::distanceTo(const Site *other, int &dx, int &dy, int &dz, bool absolutes) const
{

    dx = m_boundaries(0)->getDistanceBetween(other->x(), m_x);
    dy = m_boundaries(1)->getDistanceBetween(other->y(), m_y);
    dz = m_boundaries(2)->getDistanceBetween(other->z(), m_z);

    if (absolutes) {
        dx = std::abs(dx);
        dy = std::abs(dy);
        dz = std::abs(dz);
    }

}

uint Site::maxDistanceTo(const Site *other) const
{
    int X, Y, Z;

    this->distanceTo(other, X, Y, Z, true);

    return findLevel((uint)X, (uint)Y, (uint)Z) + 1;

}

double Site::potentialBetween(const Site *other)
{
    int X, Y, Z;

    distanceTo(other, X, Y, Z, true);

    X += Site::nNeighborsLimit();
    Y += Site::nNeighborsLimit();
    Z += Site::nNeighborsLimit();

    return DiffusionReaction::potential(X, Y, Z);
}

void Site::setDirectUpdateFlags()
{

    for (Site * neighbor : m_allNeighbors)
    {
        if (neighbor->isActive())
        {
            //This approach assumes that recursive updating of non-neighboring sites
            //WILL NOT ACTIVATE OR DEACTIVATE any sites, simply change their state,
            //and thus not interfere with any flags set here, not require flags of their own.
            for (Reaction * reaction : neighbor->siteReactions())
            {
                reaction->setDirectUpdateFlags(this);
            }

            m_affectedSites.insert(neighbor);

        }
    }

}

bool Site::hasNeighboring(int state, int range) const
{

    Site * nextNeighbor;

    for (int i = -range; i <= range; ++i)
    {
        for (int j = -range; j <= range; ++j)
        {
            for (int k = -range; k <= range; ++k)
            {

                nextNeighbor = m_neighborHood[i + Site::nNeighborsLimit()]
                        [j + Site::nNeighborsLimit()]
                        [k + Site::nNeighborsLimit()];

                if (nextNeighbor == NULL)
                {
                    continue;
                }

                else if (nextNeighbor == this)
                {
                    continue;
                }

                else if (nextNeighbor->particleState() == state ||
                         nextNeighbor->particleState() == ParticleStates::equalAs(state))
                {
                    return true;
                }

            }
        }
    }

    return false;

}


uint Site::countNeighboring(int state, int range) const
{

    Site * nextNeighbor;
    uint count = 0;

    for (int i = -range; i <= range; ++i)
    {
        for (int j = -range; j <= range; ++j)
        {
            for (int k = -range; k <= range; ++k)
            {

                nextNeighbor = m_neighborHood[i + Site::nNeighborsLimit()]
                        [j + Site::nNeighborsLimit()]
                        [k + Site::nNeighborsLimit()];

                if (nextNeighbor == NULL)
                {
                    continue;
                }

                else if (nextNeighbor == this)
                {
                    continue;
                }

                else if (nextNeighbor->particleState() == state ||
                         nextNeighbor->particleState() == ParticleStates::equalAs(state))
                {
                    count++;
                }

            }
        }
    }

    return count;

}

void Site::activate()
{

    KMCDebugger_AssertBool(!m_active, "activating active site", info());
    KMCDebugger_AssertBool(!isCrystal(), "Activating a crystal. (should always be active)", info());


    m_active = true;
    m_affectedSites.insert(this);

    informNeighborhoodOnChange(+1);

    if (isSurface())
    {
        setParticleState(ParticleStates::crystal);
    }
    else
    {
        KMCDebugger_PushImplication(this, "activation");
    }

    setDirectUpdateFlags();

    for (Reaction * reaction : m_siteReactions)
    {
        reaction->setDirectUpdateFlags(this);
    }

    KMCDebugger_MarkPartialStep("ACTIVATION COMPLETE");

    m_totalActiveSites++;


}

void Site::deactivate()
{

    KMCDebugger_AssertBool(m_active, "deactivating deactive site. ", info());
    KMCDebugger_AssertBool(!isSurface(), "deactivating a surface. (should always be deactive)", info());


    m_active = false;


    informNeighborhoodOnChange(-1);

    //if we deactivate a crystal site, we have to potentially
    //reduce the surface by adding more sites as solution sites.
    //Site will change only if it is not surrounded by any crystals.
    if (isCrystal())
    {
        KMCDebugger_Assert(particleState(), !=, ParticleStates::surface);
        KMCDebugger_Assert(particleState(), !=, ParticleStates::solution);

        setParticleState(ParticleStates::surface);
    }


    else if (qualifiesAsSurface())
    {

        KMCDebugger_Assert(particleState(), ==, ParticleStates::solution);

        m_particleState = ParticleStates::surface;

        KMCDebugger_MarkPre("surfaceSolution");
        KMCDebugger_PushImplication(this, "surface");
    }

    else
    {
        KMCDebugger_PushImplication(this, "deactivation");
    }


    setDirectUpdateFlags();


    m_activeReactions.clear();

    KMCDebugger_MarkPartialStep("DEACTIVATION COMPLETE");

    m_totalActiveSites--;

}


void Site::introduceNeighborhood()
{

    uint xTrans, yTrans, zTrans;

    Site * neighbor;

    KMCDebugger_Assert(m_nNeighborsLimit, !=, 0, "Neighborlimit must be greater than zero.", info());


    m_nNeighbors.zeros(m_nNeighborsLimit);

    Boundary::setupCurrentBoundaries(x(), y(), z());

    m_neighborHood = new Site***[m_neighborhoodLength];

    for (uint i = 0; i < m_neighborhoodLength; ++i)
    {

        xTrans = Boundary::currentBoundaries(0)->transformCoordinate((int)m_x + m_originTransformVector(i));

        m_neighborHood[i] = new Site**[m_neighborhoodLength];

        for (uint j = 0; j < m_neighborhoodLength; ++j)
        {

            yTrans = Boundary::currentBoundaries(1)->transformCoordinate((int)m_y + m_originTransformVector(j));

            m_neighborHood[i][j] = new Site*[m_neighborhoodLength];

            for (uint k = 0; k < m_neighborhoodLength; ++k)
            {

                zTrans = Boundary::currentBoundaries(2)->transformCoordinate((int)m_z + m_originTransformVector(k));

                if (Boundary::isBlocked(xTrans) ||
                        Boundary::isBlocked(yTrans) ||
                        Boundary::isBlocked(zTrans))
                {
                    m_neighborHood[i][j][k] = NULL;
                }

                else
                {

                    neighbor = m_solver->getSite(xTrans, yTrans, zTrans);

                    m_neighborHood[i][j][k] = neighbor;

                    if (neighbor != this)
                    {

                        KMCDebugger_AssertBool(!(i == Site::nNeighborsLimit() && j == Site::nNeighborsLimit() && k == Site::nNeighborsLimit()));
                        KMCDebugger_AssertBool(!(neighbor->x() == x() && neighbor->y() == y() && neighbor->z() == z()));
                        KMCDebugger_AssertBool(!(xTrans == x() && yTrans == y() && zTrans == z()));

                        m_allNeighbors.push_back(neighbor);

                        if (neighbor->isActive())
                        {

                            uint level = m_levelMatrix(i, j, k);

                            m_nNeighbors(level)++;

                            m_nNeighborsSum++;

                            double dE = DiffusionReaction::potential(i,  j,  k);

                            m_energy += dE;

                            m_totalEnergy += dE;

                        }

                    }
                }

            }
        }
    }

}

void Site::propagateToNeighbors(int reqOldState, int newState, int range)
{

    Site *nextNeighbor;

    KMCDebugger_Assert(range, <=, (int)Site::nNeighborsLimit(), "cannot propagate information beyond neighbor limit.");

    bool acceptAnything = reqOldState == ParticleStates::any;

    for (int i = -range; i <= range; ++i)
    {
        for (int j = -range; j <= range; ++j)
        {
            for (int k = -range; k <= range; ++k)
            {

                nextNeighbor = m_neighborHood[i + Site::nNeighborsLimit()]
                        [j + Site::nNeighborsLimit()]
                        [k + Site::nNeighborsLimit()];

                if (nextNeighbor == NULL)
                {
                    continue;
                }

                else if (nextNeighbor == this)
                {
                    assert(i == j && j == k && k == 0);
                    continue;
                }

                else if (nextNeighbor->particleState() == newState)
                {
                    continue;
                }

                else if (nextNeighbor->particleState() == reqOldState || acceptAnything)
                {
                    nextNeighbor->setParticleState(newState);
                }

            }
        }
    }

}

void Site::informNeighborhoodOnChange(int change)
{

    Site *neighbor;
    uint level;
    double dE;

    for (uint i = 0; i < m_neighborhoodLength; ++i)
    {
        for (uint j = 0; j < m_neighborhoodLength; ++j)
        {
            for (uint k = 0; k < m_neighborhoodLength; ++k)
            {

                neighbor = m_neighborHood[i][j][k];

                if (neighbor == NULL)
                {
                    continue;
                }

                else if (neighbor == this) {
                    assert(i == j && j == k && k == m_nNeighborsLimit);
                    continue;
                }

                KMCDebugger_AssertBool(!((change < 0) && (neighbor->nNeighborsSum() == 0)), "Call initiated to set negative nNeighbors.", neighbor->info());

                level = m_levelMatrix(i, j, k);

                KMCDebugger_AssertBool(!((change < 0) && (neighbor->nNeighbors(level) == 0)), "Call initiated to set negative neighbor.", neighbor->info());

                neighbor->m_nNeighbors(level)+= change;
                neighbor->m_nNeighborsSum += change;


                dE = change*DiffusionReaction::potential(i,  j,  k);

                neighbor->m_energy += dE;

                m_totalEnergy += dE;

            }
        }
    }


}

void Site::queueAffectedSites()
{
    for (Site * neighbor : m_allNeighbors)
    {
        if (neighbor->isActive())
        {

            for (Reaction * reaction : neighbor->siteReactions())
            {
                reaction->addUpdateFlag(Reaction::defaultUpdateFlag);
            }

            m_affectedSites.insert(neighbor);
        }
    }
}


//! Sets the site energy to zero. Used to avoid round-off error zeros.
void Site::setZeroEnergy()
{
    KMCDebugger_Assert(m_nNeighborsSum, ==, 0, "Energy is not zero.", info());
    m_energy = 0;
}

void Site::clearNeighborhood()
{

    m_totalEnergy -= m_energy;

    m_energy = 0;


    m_allNeighbors.clear();

    m_nNeighbors.reset();

    m_nNeighborsSum = 0;


    for (uint i = 0; i < m_neighborhoodLength; ++i)
    {
        for (uint j = 0; j < m_neighborhoodLength; ++j)
        {
            for (uint k = 0; k < m_neighborhoodLength; ++k)
            {
                m_neighborHood[i][j][k] = NULL;
            }

            delete [] m_neighborHood[i][j];
        }

        delete [] m_neighborHood[i];
    }

    delete [] m_neighborHood;

}

uint Site::findLevel(uint i, uint j, uint k)
{

    uint m = i;

    if (j > i)
    {
        m =  j;
    }

    if (k > m)
    {
        m = k;
    }

    return m - 1;

}

void Site::clearAll()
{

    m_nNeighborsToCrystallize = KMCSolver::UNSET_UINT;
    m_nNeighborsLimit = KMCSolver::UNSET_UINT;
    m_neighborhoodLength = KMCSolver::UNSET_UINT;

    m_totalActiveSites = 0;
    m_totalEnergy = 0;
    m_levelMatrix.reset();
    m_originTransformVector.reset();

    clearAffectedSites();
    clearBoundaries();

    m_boundaryConfigs.clear();
    m_boundaryTypes.clear();

}

void Site::resetBoundariesTo(const umat &boundaryMatrix)
{

    m_solver->clearSiteNeighborhoods();

    clearBoundaries();


    setBoundaries(boundaryMatrix);


    m_solver->initializeSiteNeighborhoods();


    m_solver->clearAllReactions();

    Site::initializeBoundaries();

    m_solver->initializeDiffusionReactions();


}

void Site::resetNNeighborsLimitTo(const uint &nNeighborsLimit)
{

    m_solver->clearSiteNeighborhoods();

    setNNeighborsLimit(nNeighborsLimit);

    m_solver->initializeSiteNeighborhoods();

}

void Site::resetNNeighborsToCrystallizeTo(const uint &nNeighborsToCrystallize)
{
    setNNeighborsToCrystallize(nNeighborsToCrystallize);
}

void Site::clearBoundaries()
{
    for (uint i = 0; i < 3; ++i)
    {
        for (uint j = 0; j < 2; ++j) {
            delete m_boundaries(i, j);
        }
    }

    m_boundaries.clear();
}

void Site::clearAffectedSites()
{
    m_affectedSites.clear();
}

void Site::setZeroTotalEnergy()
{
    m_totalEnergy = 0;
}


const string Site::info(int xr, int yr, int zr, string desc) const
{

    stringstream s_full;

    s_full << str();
    s_full << "[" << NX() << " x " << NY() << " x " << NZ() << "] * ";

    if (m_active)
    {
        s_full << "Active";
    }

    else
    {
        s_full << "Deactive";
    }

    s_full << " " << ParticleStates::names.at(particleState());

    s_full << " * Neighbors: ";

    for (uint n : m_nNeighbors)
    {
        s_full << n << " ";

    }

    s_full << "\n";

    uint _min = ParticleStates::surface + 1;

    ucube nN;
    nN.copy_size(m_levelMatrix);
    nN.fill(_min);

    Site * currentSite;
    for (uint i = 0; i < m_neighborhoodLength; ++i)
    {
        for (uint j = 0; j < m_neighborhoodLength; ++j)
        {
            for (uint k = 0; k < m_neighborhoodLength; ++k)
            {

                currentSite = m_neighborHood[i][j][k];


                if (currentSite == NULL)
                {
                    nN(i, j, k) = _min + 1;
                }

                else if (currentSite == this)
                {
                    nN(i, j, k) = _min + 2;
                }

                else if ((i == Site::nNeighborsLimit() + xr) && (j == Site::nNeighborsLimit() + yr) && (k == Site::nNeighborsLimit() + zr))
                {
                    nN(i, j, k) = _min + 3;
                }

                else if (currentSite->isActive())
                {
                    nN(i, j, k) = currentSite->particleState();
                }

                else if (currentSite->isSurface())
                {
                    nN(i, j, k) = ParticleStates::surface;
                }

            }

        }

    }

    umat A;
    stringstream ss;

    for (int j = m_neighborhoodLength - 1; j >= 0; --j)
    {
        for(uint i = 0; i < m_neighborhoodLength; ++i)
        {
            A = nN.slice(i).t();

            ss << " ";

            for (auto val : A.row(j).eval())
            {
                ss << val << " ";
            }

            if (i != m_neighborhoodLength - 1) ss << " | ";
        }

        ss << "\n";

    }


    string s = ss.str();

    auto searchRepl = [&s] (string _find, string _repl)
    {

        int position = s.find(_find);
        while (position != (int)string::npos)
        {
            s.replace(position, _find.size(), _repl);
            position = s.find(_find, position + 1);
        }

    };

    auto numberSearchRepl = [&s, &searchRepl] (int number, string desc)
    {
        stringstream type;
        type << number;
        for (uint i = 1; i < desc.size(); ++i) {
            type << " ";
        }
        searchRepl(type.str(), desc);
    };


    searchRepl("        ", "  ");

    numberSearchRepl(ParticleStates::crystal,      ParticleStates::shortNames.at(ParticleStates::crystal));
    numberSearchRepl(ParticleStates::fixedCrystal, ParticleStates::shortNames.at(ParticleStates::fixedCrystal));
    numberSearchRepl(ParticleStates::surface,      ParticleStates::shortNames.at(ParticleStates::surface));
    numberSearchRepl(ParticleStates::solution,     ParticleStates::shortNames.at(ParticleStates::solution));

    numberSearchRepl(_min+0, ".");
    numberSearchRepl(_min+1, "--");
    numberSearchRepl(_min+2, particleStateShortName() + "^");
    numberSearchRepl(_min+3, desc);

    s_full << s;

    string full_string = s_full.str();

    return full_string;

}

uint Site::nNeighborsSum() const
{
    KMCDebugger_Assert(sum(m_nNeighbors), ==, m_nNeighborsSum, "Should be identical.", info());
    return m_nNeighborsSum;
}

void Site::setNNeighborsLimit(const uint &nNeighborsLimit)
{

    if (nNeighborsLimit >= min(uvec({NX(), NY(), NZ()}))/2)
    {
        cerr << "Neighbor reach must be lower than half the minimum box dimension to avoid sites directly affecting themselves." << endl;
        exit(1);
    }


    m_nNeighborsLimit = nNeighborsLimit;

    m_neighborhoodLength = 2*m_nNeighborsLimit + 1;

    m_levelMatrix.set_size(m_neighborhoodLength, m_neighborhoodLength, m_neighborhoodLength);

    m_originTransformVector = linspace<ivec>(-(int)m_nNeighborsLimit, m_nNeighborsLimit, m_neighborhoodLength);

    for (uint i = 0; i < m_neighborhoodLength; ++i)
    {
        for (uint j = 0; j < m_neighborhoodLength; ++j)
        {
            for (uint k = 0; k < m_neighborhoodLength; ++k)
            {
                if (i == m_nNeighborsLimit && j == m_nNeighborsLimit && k == m_nNeighborsLimit)
                {
                    m_levelMatrix(i, j, k) = m_nNeighborsLimit + 1;
                    continue;
                }

                m_levelMatrix(i, j, k) = findLevel(std::abs(m_originTransformVector(i)),
                                                   std::abs(m_originTransformVector(j)),
                                                   std::abs(m_originTransformVector(k)));
            }
        }
    }

    DiffusionReaction::setupPotential();

}

void Site::setNNeighborsToCrystallize(const uint &nNeighborsToCrystallize)

{
    if (nNeighborsToCrystallize == 0)
    {
        cerr << "With nNeighborsToCrystallize = 0, all particles will qualify as crystals." << endl;
        exit(1);
    }
    else if (nNeighborsToCrystallize > 7)
    {
        cerr << "With nNeighboorsToCrystallize > 7, no particles except those hugging a fixed crystal will qualify as crystals." << endl;
        exit(1);
    }

    m_nNeighborsToCrystallize = nNeighborsToCrystallize;

}

void Site::setInitialBoundaries(const Setting & boundariesConfig)
{


    m_boundaryConfigs.set_size(3, 2);
    m_boundaryTypes.set_size(3, 2);

    for (uint XYZ = 0; XYZ < 3; ++XYZ)
    {
        for (uint orientation = 0; orientation < 2; ++orientation)
        {
            m_boundaryConfigs(XYZ, orientation) = &getSurfaceSetting(boundariesConfig, "configs")[XYZ][orientation];
            m_boundaryTypes(XYZ, orientation) = getSurfaceSetting(boundariesConfig, "types")[XYZ][orientation];
        }
    }

    setBoundaries(m_boundaryTypes);

}

void Site::setBoundaries(const umat &boundaryMatrix)
{

    m_boundaryTypes = boundaryMatrix;

    m_boundaries.set_size(3, 2);

    for (uint XYZ = 0; XYZ < 3; ++XYZ)
    {
        for (uint orientation = 0; orientation < 2; ++orientation)
        {

            switch (m_boundaryTypes(XYZ, orientation))
            {
            case Boundary::Periodic:
                m_boundaries(XYZ, orientation) = new Periodic(XYZ, orientation);

                break;

            case Boundary::Edge:
                m_boundaries(XYZ, orientation) = new Edge(XYZ, orientation);

                break;

            case Boundary::Surface:
                m_boundaries(XYZ, orientation) = new Surface(XYZ, orientation);

                break;

            case Boundary::ConcentrationWall:
                m_boundaries(XYZ, orientation) = new ConcentrationWall(XYZ, orientation);

                break;

            default:

                cerr << "Unknown boundary type " << m_boundaryTypes(XYZ, orientation) << endl;
                exit(1);

                break;
            }

            m_boundaries(XYZ, orientation)->loadConfig(*m_boundaryConfigs(XYZ, orientation));


        }

        if (!Boundary::isCompatible(m_boundaryTypes(XYZ, 0), m_boundaryTypes(XYZ, 1)))
        {
            cerr << "Mismatch in boundaries for " << XYZ << "'th dimension: " << m_boundaryTypes.t();
            exit(1);
        }
    }

}



const uint &Site::NX()
{
    return m_solver->NX();
}

const uint &Site::NY()
{
    return m_solver->NY();
}

const uint &Site::NZ()
{
    return m_solver->NZ();
}







KMCSolver* Site::m_solver;

uint       Site::m_nNeighborsToCrystallize = KMCSolver::UNSET_UINT;

uint       Site::m_nNeighborsLimit = KMCSolver::UNSET_UINT;

uint       Site::m_neighborhoodLength = KMCSolver::UNSET_UINT;

ucube      Site::m_levelMatrix;
ivec       Site::m_originTransformVector;

uint       Site::m_totalActiveSites = 0;

double     Site::m_totalEnergy = 0;

set<Site*> Site::m_affectedSites;

field<Boundary*> Site::m_boundaries;

field<const Setting*> Site::m_boundaryConfigs;

umat Site::m_boundaryTypes;

const vector<string> ParticleStates::names = {"crystal", "fixedcrystal", "solution", "surface"};
const vector<string> ParticleStates::shortNames = {"C", "F", "P", "S"};


ostream & operator << (ostream& os, const Site& ss)
{
    os << ss.str();
    return os;
}


int ParticleStates::equalAs(int state)
{
    switch (state) {
    case ParticleStates::crystal:
        return ParticleStates::fixedCrystal;
        break;
    default:
        return state;
        break;
    }
}



