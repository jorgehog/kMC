#include "site.h"
#include "kmcsolver.h"
#include "reactions/reaction.h"
#include "reactions/diffusion/diffusionreaction.h"

#include "boundary/periodic/periodic.h"
#include "boundary/edge/edge.h"
#include "boundary/wall/wall.h"
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
    m_energy(0),
    m_particleState(ParticleStates::solution)
{

}


Site::~Site()
{
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


    for (Reaction* reaction : m_siteReactions)
    {
        delete reaction;
    }

    m_activeReactions.clear();

    m_siteReactions.clear();

    m_allNeighbors.clear();

    m_nNeighbors.reset();

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


void Site::setParticleState(int state)
{

    //If we try to propagate a surface onto an existing active site
    //The site is crystallized, which propagates the surface furhter.
    switch (state)
    {
    case ParticleStates::surface:

        switch (m_particleState)
        {

        //solution->surface
        case ParticleStates::solution:

            //if a particle is present, we crystallize it if it fulfills the crystallizaton criteria.
            if (m_active)
            {
                if (shouldCrystallize())
                {
                    m_particleState  = ParticleStates::crystal;
                    KMCDebugger_PushImplication(this, "solution", "crystal");
                    propagateToNeighbors(ParticleStates::solution, ParticleStates::surface, DiffusionReaction::separation());
                }
            }

            else
            {
                m_particleState = ParticleStates::surface;
                queueAffectedSites();
                KMCDebugger_PushImplication(this, "solution", "surface");
            }


            break;

            //Crystal -> surface
        case (ParticleStates::crystal):

            if (hasNeighboring(ParticleStates::crystal, DiffusionReaction::separation()))
            {
                m_particleState = ParticleStates::surface;
                KMCDebugger_PushImplication(this, "crystal", "surface");
            }

            else
            {
                m_particleState = ParticleStates::solution;
                KMCDebugger_PushImplication(this, "crystal", "solution");
            }

            propagateToNeighbors(ParticleStates::surface, ParticleStates::solution, DiffusionReaction::separation());
            queueAffectedSites();

            break;

            //surface -> surface
        case ParticleStates::surface:
            //Nothing to do here.
            KMCDebugger_PushImplication(this, "surface", "surface");
            break;

        default:
            cout << "invalid transition "
                 << ParticleStates::names.at(m_particleState)
                 << "->"
                 << ParticleStates::names.at(state) << endl;
            exit(1);
            break;
        }

        break;

    case ParticleStates::crystal:

        switch (m_particleState)
        {

        //surface -> crystal
        case ParticleStates::surface:

            //No need to test if it has neigh crystals because
            //diffusion reactions deactivates old spot before activating new spot.
            //Which will remove access surface.
            m_particleState  = ParticleStates::crystal;
            KMCDebugger_PushImplication(this, "surface", "crystal");
            propagateToNeighbors(ParticleStates::solution, ParticleStates::surface, DiffusionReaction::separation());

            break;

            //fixedCrystal -> crystal
        case ParticleStates::fixedCrystal:

            cout << "Warning: Making fixed crystal seed into standard crystal." << endl;
            m_particleState = ParticleStates::crystal;
            KMCDebugger_PushImplication(this, "fixedcrystal", "crystal");
            break;

        default:
            cout << "invalid transition "
                 << ParticleStates::names.at(m_particleState)
                 << "->"
                 << ParticleStates::names.at(state) << endl;
            exit(1);
            break;
        }

        break;


    case ParticleStates::solution:

        switch (m_particleState)
        {

        //surface -> solution
        case ParticleStates::surface:

            if (!(hasNeighboring(ParticleStates::crystal, DiffusionReaction::separation()) || m_isFixedCrystalSeed))
            {
                m_particleState = ParticleStates::solution;
                KMCDebugger_PushImplication(this, "surface", "solution");
                queueAffectedSites();
            }
#ifndef KMC_NO_DEBUG
            else
            {
                KMCDebugger_PushImplication(this, "surface", "surface");
            }
#endif

            break;


        default:
            cout << "invalid transition "
                 << ParticleStates::names.at(m_particleState)
                 << "->"
                 << ParticleStates::names.at(state) << endl;
            exit(1);
            break;
        }
        break;



    default:
        cout << "invalid transition "
             << ParticleStates::names.at(m_particleState)
             << "->"
             << ParticleStates::names.at(state) << endl;
        exit(1);
        break;
    }

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

bool Site::shouldCrystallize()
{

    KMCDebugger_AssertBool((isActive() == true) && (particleState() == ParticleStates::solution), "only active solution particles should crystallize.", info());

    //Dummy
    return true;

}


void Site::loadConfig(const Setting &setting)
{

    m_boundaries.set_size(3, 2);

    const Setting & boundariesConfig = getSurfaceSetting(setting, "Boundaries");

    ivec boundaryTypes(2);

    for (uint XYZ = 0; XYZ < 3; ++XYZ)
    {
        for (uint orientation = 0; orientation < 2; ++orientation)
        {

            boundaryTypes(orientation) = getSurfaceSetting(boundariesConfig, "types")[XYZ][orientation];

            switch (boundaryTypes(orientation))
            {
            case Boundary::Periodic:
                m_boundaries(XYZ, orientation) = new Periodic(XYZ, orientation);

                break;

            case Boundary::Edge:
                m_boundaries(XYZ, orientation) = new Edge(XYZ, orientation);

                break;

            case Boundary::Wall:
                m_boundaries(XYZ, orientation) = new Wall(XYZ, orientation);

                break;

            case Boundary::ConsentrationWall:
                m_boundaries(XYZ, orientation) = new ConcentrationWall(XYZ, orientation);

                break;

            default:

                cerr << "Unknown boundary type " << boundaryTypes(orientation) << endl;
                exit(1);

                break;
            }

            m_boundaries(XYZ, orientation)->loadConfig(getSurfaceSetting(boundariesConfig, "configs")[XYZ][orientation]);

        }

        if (!Boundary::isCompatible(boundaryTypes(0), boundaryTypes(1)))
        {
            cerr << "Mismatch in boundaries for " << XYZ << "'th dimension: " << boundaryTypes.t();
            exit(1);
        }


    }


    const uint  &limit = getSurfaceSetting<uint>(setting, "nNeighborsLimit");

    if (limit >= min(uvec({NX, NY, NZ}))/2)
    {
        cerr << "Neighbor reach must be lower than half the minimum box dimension to avoid sites directly affecting themselves." << endl;
        exit(1);
    }

    m_nNeighborsLimit = limit;
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

}

void Site::initializeBoundaries()
{
    for (uint i = 0; i < 3; ++i) {
        for (uint j = 0; j < 2; ++j) {
            m_boundaries(i, j)->initialize();
        }
    }
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


void Site::spawnAsFixedCrystal()
{
    m_particleState = ParticleStates::surface;
    m_isFixedCrystalSeed = true;

    for (Reaction * reaction : m_siteReactions)
    {
        delete reaction;
    }

    m_siteReactions.clear();

    activate();

    m_particleState = ParticleStates::fixedCrystal;

}


void Site::calculateRates()
{
    for (Reaction* reaction : m_activeReactions)
    {
        reaction->selectTriumphingUpdateFlag();
        reaction->calcRate();
    }
}


void Site::setMainSolver(KMCSolver *solver)
{

    mainSolver = solver;

    NX = solver->getNX();
    NY = solver->getNY();
    NZ = solver->getNZ();

}


void Site::distanceTo(const Site *other, int &dx, int &dy, int &dz, bool absolutes) const
{

    dx = m_boundaries(0)->getDistanceBetween(other->x(), m_x);
    dy = m_boundaries(1)->getDistanceBetween(other->y(), m_y);
    dz = m_boundaries(2)->getDistanceBetween(other->z(), m_z);

    if (absolutes) {
        dx = abs(dx);
        dy = abs(dy);
        dz = abs(dz);
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

void Site::activate()
{

    KMCDebugger_AssertBool(!m_active, "activating deactivated site", info());
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
        KMCDebugger_PushImplication(this, "deactiveSolution", "activeSolution");
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
        setParticleState(ParticleStates::surface);
    }

    else
    {
        KMCDebugger_PushImplication(this, "activeSolution", "deactiveSolution");
    }

    setDirectUpdateFlags();


    m_activeReactions.clear();

    KMCDebugger_MarkPartialStep("DEACTIVATION COMPLETE");

    m_totalActiveSites--;

}


void Site::introduceNeighborhood()
{

    assert(m_nNeighborsLimit != 0);

    uint xTrans, yTrans, zTrans;
    uvec3 loc;

    Boundary::setupLocations(m_x, m_y, m_z, loc);

    const Boundary * xBoundary = m_boundaries(0, loc(0));
    const Boundary * yBoundary = m_boundaries(1, loc(1));
    const Boundary * zBoundary = m_boundaries(2, loc(2));


    m_nNeighbors.set_size(m_nNeighborsLimit);
    m_nNeighbors.zeros();

    m_neighborHood = new Site***[m_neighborhoodLength];


    for (uint i = 0; i < m_neighborhoodLength; ++i)
    {

        xTrans = xBoundary->transformCoordinate((int)m_x + m_originTransformVector(i));

        m_neighborHood[i] = new Site**[m_neighborhoodLength];

        for (uint j = 0; j < m_neighborhoodLength; ++j)
        {

            yTrans = yBoundary->transformCoordinate((int)m_y + m_originTransformVector(j));

            m_neighborHood[i][j] = new Site*[m_neighborhoodLength];

            for (uint k = 0; k < m_neighborhoodLength; ++k)
            {

                zTrans = zBoundary->transformCoordinate((int)m_z + m_originTransformVector(k));

                if (Boundary::isBlocked(xTrans) ||
                        Boundary::isBlocked(yTrans) ||
                        Boundary::isBlocked(zTrans))
                {
                    m_neighborHood[i][j][k] = NULL;
                }

                else
                {

                    m_neighborHood[i][j][k] = mainSolver->getSite(xTrans, yTrans, zTrans);

                    if (m_neighborHood[i][j][k] != this)
                    {
                        m_allNeighbors.push_back(m_neighborHood[i][j][k]);
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

                KMCDebugger_AssertBool(!(newState == ParticleStates::solution && nextNeighbor->particleState() == ParticleStates::solution),
                                       "Solution asked to go to solution.", info());
                if (nextNeighbor->particleState() == reqOldState)
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

void Site::resetAll()
{

    m_totalActiveSites = 0;
    m_totalEnergy = 0;
    m_levelMatrix.reset();
    m_originTransformVector.reset();
    m_affectedSites.clear();

    for (uint i = 0; i < 3; ++i)
    {
        for (uint j = 0; j < 2; ++j) {
            delete m_boundaries(i, j);
        }
    }

    m_boundaries.clear();

}


const string Site::info(int xr, int yr, int zr, string desc) const
{

    stringstream s_full;

    s_full << str();
    s_full << "[" << NX << " x " << NY << " x " << NZ << "] * ";

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
    numberSearchRepl(_min+1, " ");
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


KMCSolver* Site::mainSolver;

uint       Site::NX;
uint       Site::NY;
uint       Site::NZ;

uint       Site::m_nNeighborsLimit;

uint       Site::m_neighborhoodLength;

ucube      Site::m_levelMatrix;
ivec       Site::m_originTransformVector;

uint       Site::m_totalActiveSites = 0;

double     Site::m_totalEnergy = 0;

set<Site*> Site::m_affectedSites;

field<Boundary*> Site::m_boundaries;

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
