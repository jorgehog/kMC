#include "site.h"
#include "reactions/reaction.h"
#include "reactions/diffusion/diffusionreaction.h"
#include "kmcsolver.h"

#include "debugger/kmcdebugger.h"


Site::Site(uint _x, uint _y, uint _z) :
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
    for (uint i = 0; i < neighborhoodLength(); ++i)
    {
        for (uint j = 0; j < neighborhoodLength(); ++j)
        {
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

    m_energy = 0;

}




void Site::updateAffectedSites()
{

    for (Site* site : affectedSites)
    {
        site->updateReactions();
        site->calculateRates();
    }

    affectedSites.clear();

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

            //if a particle is present, we crystallize it immidiately.
            if (m_active)
            {
                m_particleState  = ParticleStates::crystal;
                KMCDebugger_PushImplication(this, "solution", "crystal");
                propagateToNeighbors(ParticleStates::solution, ParticleStates::surface);
            }

            else
            {
                m_particleState = ParticleStates::surface;
                KMCDebugger_PushImplication(this, "solution", "surface");
            }

            queueAffectedSites();

            break;

            //Crystal -> surface
        case ParticleStates::crystal:

            if (m_isFixedCrystalSeed)
            {
                m_particleState = ParticleStates::surface;
                KMCDebugger_PushImplication(this, "fixedcrystal", "surface");
            }

            else if (hasNeighboring(ParticleStates::crystal))
            {
                m_particleState = ParticleStates::surface;
                KMCDebugger_PushImplication(this, "crystal", "surface");
            }

            else
            {
                m_particleState = ParticleStates::solution;
                KMCDebugger_PushImplication(this, "crystal", "solution");
            }

            propagateToNeighbors(ParticleStates::surface, ParticleStates::solution);
            queueAffectedSites();

            break;

            //surface -> surface
        case ParticleStates::surface:
            //Nothing to do here.
            KMCDebugger_PushImplication(this, "surface", "surface");
            break;

        default:
            cout << "invalid transition"
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
            propagateToNeighbors(ParticleStates::solution, ParticleStates::surface);
            queueAffectedSites();

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

            if (!(hasNeighboring(ParticleStates::crystal) || m_isFixedCrystalSeed))
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
        if (!r->allowedAtSite())
        {
            return false;
        }
    }


    return true;

}


void Site::loadConfig(const Setting &setting)
{

    const uint  &limit = getSurfaceSetting<uint>(setting, "nNeighborsLimit");

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


void Site::addReaction(Reaction *reaction)
{
    reaction->setSite(this);

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
        if (reaction->isNotBlocked())
        {
            m_activeReactions.push_back(reaction);
        }
    }

}


void Site::spawnAsFixedCrystal()
{
    m_particleState = ParticleStates::surface;
    m_isFixedCrystalSeed = true;
    activate();
}


void Site::calculateRates()
{
    for (Reaction* reaction : m_activeReactions)
    {
        reaction->getTriumphingUpdateFlag();
        reaction->calcRate();
    }
}


void Site::setSolverPtr(KMCSolver *solver)
{

    NX = solver->getNX();
    NY = solver->getNY();
    NZ = solver->getNZ();

    mainSolver = solver;

}


void Site::distanceTo(const Site *other, int &dx, int &dy, int &dz, bool absolutes) const
{

    dx = (other->x() + NX - m_x)%NX;
    dy = (other->y() + NY - m_y)%NY;
    dz = (other->z() + NZ - m_z)%NZ;


    if ((uint)abs(dx) > NX/2) {
        dx = -(int)(NX - dx);
    }

    if ((uint)abs(dy) > NY/2) {
        dy = -(int)(NY - dy);
    }

    if ((uint)abs(dz) > NZ/2) {
        dz = -(int)(NZ - dz);
    }


    if (absolutes) {
        dx = abs(dx);
        dy = abs(dy);
        dz = abs(dz);
    }

}

uint Site::maxDistanceTo(const Site *other)
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

bool Site::hasNeighboring(int state) const
{

    Site *nextNeighbor;

    for (uint i = 0; i < 3; ++i)
    {
        for (uint j = 0; j < 3; ++j)
        {
            for (uint k = 0; k < 3; ++k)
            {
                if (i == 1 && j == 1 && k == 1)
                {
                    continue;
                }

                nextNeighbor = m_neighborHood[i + Site::nNeighborsLimit() - 1]
                        [j + Site::nNeighborsLimit() - 1]
                        [k + Site::nNeighborsLimit() - 1];

                if (nextNeighbor->particleState() == state)
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

#ifndef KMC_NO_DEBUG

    if (m_active == true)
    {
        cout << "Activating active site. " << endl;
        info();
        exit(1);
    }

    else if (isCrystal())
    {
        cout << "Activating a crystal. (should always be active)";
        info();
        exit(1);
    }

#endif

    m_active = true;

    affectedSites.insert(this);

    if (isSurface())
    {
        setParticleState(ParticleStates::crystal);
    }
#ifndef KMC_NO_DEBUG
    else
    {
        KMCDebugger_PushImplication(this, "deactiveSolution", "activeSolution");
    }
#endif

    informNeighborhoodOnChange(+1);

    m_totalActiveSites++;

    KMCDebugger_MarkPartialStep("ACTIVATE_DONE");

}

void Site::deactivate()
{

#ifndef KMC_NO_DEBUG

    if (m_active == false)
    {
        cout << "deactivating deactive site. " << endl;
        info();
        exit(1);
    }
    else if (isSurface())
    {
        cout << "deactivating a surface. (should always be deactive)";
        info();
        exit(1);
    }
#endif

    m_active = false;

    //if we deactivate a crystal site, we have to potentially
    //reduce the surface by adding more sites as solution sites.
    //Site will change only if it is not surrounded by any crystals.
    if (isCrystal())
    {
        setParticleState(ParticleStates::surface);
    }
#ifndef KMC_NO_DEBUG
    else
    {
        KMCDebugger_PushImplication(this, "activeSolution", "deactiveSolution");
    }
#endif

    informNeighborhoodOnChange(-1);

    m_totalActiveSites--;

    KMCDebugger_MarkPartialStep("DEACTIVATE_DONE");

}


void Site::introduceNeighborhood()
{

    assert(m_nNeighborsLimit != 0);

    uint xTrans, yTrans, zTrans;

    m_nNeighbors.set_size(m_nNeighborsLimit);
    m_nNeighbors.zeros();

    m_neighborHood = new Site***[m_neighborhoodLength];

    for (uint i = 0; i < m_neighborhoodLength; ++i)
    {
        xTrans = (m_x + m_originTransformVector(i) + NX)%NX;

        m_neighborHood[i] = new Site**[m_neighborhoodLength];

        for (uint j = 0; j < m_neighborhoodLength; ++j)
        {
            yTrans = (m_y + m_originTransformVector(j) + NY)%NY;

            m_neighborHood[i][j] = new Site*[m_neighborhoodLength];

            for (uint k = 0; k < m_neighborhoodLength; ++k)
            {
                zTrans = (m_z + m_originTransformVector(k) + NZ)%NZ;

                m_neighborHood[i][j][k] = mainSolver->getSite(xTrans, yTrans, zTrans);

                if (m_neighborHood[i][j][k] != this)
                {
                    m_allNeighbors.push_back(m_neighborHood[i][j][k]);
                }
            }
        }
    }

}

void Site::propagateToNeighbors(int reqOldState, int newState)
{

    Site *nextNeighbor;

    for (uint i = 0; i < 3; ++i)
    {
        for (uint j = 0; j < 3; ++j)
        {
            for (uint k = 0; k < 3; ++k)
            {

                nextNeighbor = m_neighborHood[i + Site::nNeighborsLimit() - 1]
                        [j + Site::nNeighborsLimit() - 1]
                        [k + Site::nNeighborsLimit() - 1];

                if (nextNeighbor == this)
                {
                    assert(i == j && j == k && k == 1);
                    continue;
                }

                assert(!(newState == ParticleStates::solution && nextNeighbor->particleState() == ParticleStates::solution));
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

                if (neighbor == this) {
                    assert(i == j && j == k && k == m_nNeighborsLimit);
                    continue;
                }

                level = m_levelMatrix(i, j, k);
                neighbor->m_nNeighbors(level)+=change;


                dE = change*DiffusionReaction::potential(i,  j,  k);

                neighbor->m_energy += dE;

                m_totalEnergy += dE;


            }
        }
    }

    queueAffectedSites();

}

void Site::queueAffectedSites()
{
    uint C = 0;
    for (uint i = 0; i < m_neighborhoodLength; ++i)
    {
        for (uint j = 0; j < m_neighborhoodLength; ++j)
        {
            for (uint k = 0; k < m_neighborhoodLength; ++k)
            {

                //This approach assumes that recursive updating of non-neighboring sites
                //WILL NOT ACTIVATE OR DEACTIVATE any sites, simply change their state,
                //and thus not interfere with any flags set here, not require flags of their own.
                for (Reaction * reaction : m_neighborHood[i][j][k]->siteReactions())
                {
                    reaction->setUpdateFlags(this, m_levelMatrix(i, j, k));
                    C++;
                }

            }
        }
    }

    assert(C == 125*26);

    affectedSites.insert(m_allNeighbors.begin(), m_allNeighbors.end());
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


const string Site::info(int xr, int yr, int zr, string desc) const
{

    stringstream s_full;

    s_full << "Site   " << m_x << " " << m_y << " " << m_z << endl;
    s_full << "in Box " << NX << " " << NY << " " << NZ << endl;
    s_full << "nNeighbors : " << m_nNeighbors.t();
    s_full << "type: " << ParticleStates::names.at(m_particleState) << endl;

    if (m_active)
    {
        s_full << "ACTIVE";
    }

    else
    {
        s_full << "DEACTIVE";
    }

    s_full << endl;

    ucube nN;
    nN.copy_size(m_levelMatrix);
    nN.fill(7);

    Site * currentSite;
    for (uint i = 0; i < m_neighborhoodLength; ++i)
    {
        for (uint j = 0; j < m_neighborhoodLength; ++j)
        {
            for (uint k = 0; k < m_neighborhoodLength; ++k)
            {

                currentSite = m_neighborHood[i][j][k];

                if (currentSite->isFixedCrystalSeed())
                {
                    assert((currentSite->isCrystal() && currentSite->isActive()) || (currentSite->isSurface() && !currentSite->isActive()));
                }

                if (currentSite == this)
                {
                    nN(i, j, k) = 9;
                }

                else if ((i == Site::nNeighborsLimit() + xr) && (j == Site::nNeighborsLimit() + yr) && (k == Site::nNeighborsLimit() + zr))
                {
                    nN(i, j, k) = 8;
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
    for(int i = nN.n_slices - 1; i >= 0; --i)
    {
        A = nN.slice(i).t();

        for (int j = A.n_rows - 1; j >= 0; --j)
        {
            ss << A.row(j);
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

    auto typeSearchRepl = [&s, &searchRepl] (int pType)
    {
        stringstream type;
        type << pType;
        searchRepl(type.str(), ParticleStates::shortNames.at(pType));
    };


    searchRepl("        ", "  ");

    typeSearchRepl(ParticleStates::crystal);
    typeSearchRepl(ParticleStates::surface);
    typeSearchRepl(ParticleStates::solution);

    searchRepl("9 ", particleStateShortName() + "^");
    searchRepl("8", desc);
    searchRepl("7", ".");

    s_full << s;

    return s_full.str();

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

set<Site*> Site::affectedSites;


const vector<string> ParticleStates::names = {"crystal", "solution", "surface"};
const vector<string> ParticleStates::shortNames = {"C", "P", "S"};


ostream & operator << (ostream& os, const Site& ss)
{
    os << ss.str();
    return os;
}
