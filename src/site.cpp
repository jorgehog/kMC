#include "site.h"
#include "kmcsolver.h"
#include "reactions/reaction.h"
#include "reactions/diffusion/diffusionreaction.h"

#include "boundary/periodic/periodic.h"
#include "boundary/edge/edge.h"
#include "boundary/concentrationwall/concentrationwall.h"

#include "soluteparticle.h"

#include "debugger/debugger.h"

using namespace kMC;


Site::Site(uint _x, uint _y, uint _z) :
    m_associatedParticle(NULL),
    m_x(_x),
    m_y(_y),
    m_z(_z)
{
    refCounter++;
}


Site::~Site()
{

    clearNeighborhood();

    refCounter--;

}


void Site::loadConfig(const Setting &setting)
{

    setInitialNNeighborsLimit(getSurfaceSetting<uint>(setting, "nNeighborsLimit"));

    const Setting & boundariesConfig = getSurfaceSetting(setting, "Boundaries");

    umat boundaryTypes(3, 2);

    for (uint XYZ = 0; XYZ < 3; ++XYZ)
    {
        for (uint orientation = 0; orientation < 2; ++orientation)
        {
            boundaryTypes(XYZ, orientation) = getSurfaceSetting(boundariesConfig, "types")[XYZ][orientation];
        }
    }

    setInitialBoundaries(boundaryTypes);

}

void Site::initializeBoundaries()
{
    KMCDebugger_Assert(Site::_refCount(), !=, 0, "Sites needs to be enabled to initialize boundaries.");
    KMCDebugger_AssertBool(!Site::boundariesIsInitialized(), "Boundaries needs to be finalized before they can be initialized.");

    KMCDebugger_SetEnabledTo(false);

    for (uint i = 0; i < 3; ++i) {
        for (uint j = 0; j < 2; ++j) {
            m_boundaries(i, j)->initialize();
            m_boundaries(i, j)->setAsInitialized();
        }
    }

    KMCDebugger_ResetEnabled();
}

void Site::updateBoundaries()
{
    for (uint i = 0; i < 3; ++i) {
        for (uint j = 0; j < 2; ++j) {
            m_boundaries(i, j)->update();
        }
    }
}

bool Site::boundariesIsInitialized()
{

    for (uint i = 0; i < 3; ++i) {
        for (uint j = 0; j < 2; ++j) {
            if (!m_boundaries(i, j)->initialized())
            {
                return false;
            }
        }
    }

    return true;
}


uint Site::maxNeighbors()
{
    return std::pow(m_neighborhoodLength, 3) - 1;
}



void Site::forEachNeighborDo(function<void (Site *)> applyFunction) const
{

    Site * neighbor;

    for (uint i = 0; i < m_neighborhoodLength; ++i)
    {
        for (uint j = 0; j < m_neighborhoodLength; ++j)
        {
            for (uint k = 0; k < m_neighborhoodLength; ++k)
            {

                neighbor = neighborhood(i, j, k);

                if (neighbor == NULL)
                {
                    continue;
                }

                else if (neighbor == this) {
                    assert(i == j && j == k && k == m_nNeighborsLimit);
                    continue;
                }


                applyFunction(neighbor);


            }
        }
    }
}

void Site::forEachNeighborDo_sendIndices(function<void (Site *, uint, uint, uint)> applyFunction) const
{

    Site * neighbor;

    for (uint i = 0; i < m_neighborhoodLength; ++i)
    {
        for (uint j = 0; j < m_neighborhoodLength; ++j)
        {
            for (uint k = 0; k < m_neighborhoodLength; ++k)
            {

                neighbor = neighborhood(i, j, k);

                if (neighbor == NULL)
                {
                    continue;
                }

                else if (neighbor == this) {
                    assert(i == j && j == k && k == m_nNeighborsLimit);
                    continue;
                }


                applyFunction(neighbor, i, j, k);


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

    return getLevel((uint)X, (uint)Y, (uint)Z) + 1;

}



bool Site::hasNeighboring(const int state) const
{

    Site * nextNeighbor;

    for (int i = -1; i <= 1; ++i)
    {
        for (int j = -1; j <= 1; ++j)
        {
            for (int k = -1; k <= 1; ++k)
            {

                nextNeighbor = neighborhood(i + Site::nNeighborsLimit(),
                                            j + Site::nNeighborsLimit(),
                                            k + Site::nNeighborsLimit());

                if (nextNeighbor == NULL)
                {
                    continue;
                }

                else if (nextNeighbor == this)
                {
                    continue;
                }

                else if (nextNeighbor->isActive())
                {
                    if (m_associatedParticle->particleState() == state)
                    {
                        return true;
                    }
                }

            }
        }
    }

    return false;

}


uint Site::countNeighboring(int state) const
{

    Site * nextNeighbor;

    uint count = 0;

    for (int i = -1; i <= 1; ++i)
    {
        for (int j = -1; j <= 1; ++j)
        {
            for (int k = -1; k <= 1; ++k)
            {

                nextNeighbor = neighborhood(i + Site::nNeighborsLimit(),
                                            j + Site::nNeighborsLimit(),
                                            k + Site::nNeighborsLimit());

                if (nextNeighbor == NULL)
                {
                    continue;
                }

                else if (nextNeighbor == this)
                {
                    continue;
                }

                else if (nextNeighbor->isActive())
                {
                    if (m_associatedParticle->particleState() == state)
                    {
                        count++;
                    }
                }

            }
        }
    }

    return count;

}


void Site::introduceNeighborhood()
{

    KMCDebugger_Assert(m_nNeighborsLimit, !=, 0, "Neighborlimit must be greater than zero.", info());
    KMCDebugger_Assert(m_nNeighborsLimit, !=, KMCSolver::UNSET_UINT, "Neighborlimit is not set.", str());


    int n_x, n_y, n_z;

    m_neighborhood = new Site***[m_neighborhoodLength];

    for (uint i = 0; i < m_neighborhoodLength; ++i)
    {
        n_x = (int)m_x + m_originTransformVector(i);

        m_neighborhood[i] = new Site**[m_neighborhoodLength];

        for (uint j = 0; j < m_neighborhoodLength; ++j)
        {
            n_y = (int)m_y + m_originTransformVector(j);

            m_neighborhood[i][j] = new Site*[m_neighborhoodLength];

            for (uint k = 0; k < m_neighborhoodLength; ++k)
            {
                n_z = (int)m_z + m_originTransformVector(k);

                m_neighborhood[i][j][k] = m_solver->getSite(n_x, n_y, n_z);

            }
        }
    }

}



void Site::clearNeighborhood()
{

    for (uint i = 0; i < m_neighborhoodLength; ++i)
    {
        for (uint j = 0; j < m_neighborhoodLength; ++j)
        {
            for (uint k = 0; k < m_neighborhoodLength; ++k)
            {
                m_neighborhood[i][j][k] = NULL;
            }

            delete [] m_neighborhood[i][j];
        }

        delete [] m_neighborhood[i];
    }

    delete [] m_neighborhood;

}

uint Site::getLevel(uint i, uint j, uint k)
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

umat Site::getCurrentCrystalBoxTopology()
{

    umat boxTop(3, 2);

    boxTop.col(0) = ucolvec({NX(), NY(), NZ()});
    boxTop.col(1).zeros();

    for (SoluteParticle *particle : m_solver->particles())
    {
        if (particle->isCrystal())
        {

            for (uint xi = 0; xi < 3; ++xi)
            {

                if (particle->r(xi) < boxTop(xi, 0))
                {
                    boxTop(xi, 0) = particle->r(xi);
                }

                if (particle->r(xi) > boxTop(xi, 1))
                {
                    boxTop(xi, 1) = particle->r(xi);
                }

            }

        }
    }

    return boxTop;
}

void Site::clearAll()
{

    m_nNeighborsLimit = KMCSolver::UNSET_UINT;
    m_neighborhoodLength = KMCSolver::UNSET_UINT;

    m_levelMatrix.reset();
    m_originTransformVector.reset();

    clearBoundaries();

    m_boundaryConfigs.clear();
    m_boundaryTypes.clear();

}

void Site::resetBoundariesTo(const umat &boundaryMatrix)
{

    KMCDebugger_Assert(refCounter, ==, 0, "Sites must be cleared before the neighborhood length can change.");

    finalizeBoundaries();

    clearBoundaries();

    setInitialBoundaries(boundaryMatrix);

}

void Site::resetBoundariesTo(const int boundaryType)
{
    resetBoundariesTo(umat::fixed<3, 2>(fill::zeros) + boundaryType);
}

void Site::resetNNeighborsLimitTo(const uint &nNeighborsLimit, bool check)
{

    KMCDebugger_Assert(refCounter, ==, 0, "Sites must be cleared before the neighborhood length can change.");

    setInitialNNeighborsLimit(nNeighborsLimit, check);

}


void Site::clearBoundaries()
{

    for (uint i = 0; i < 3; ++i)
    {
        for (uint j = 0; j < 2; ++j)
        {
            delete m_boundaries(i, j);
        }
    }

    m_boundaries.clear();
}

void Site::finalizeBoundaries()
{
    for (uint i = 0; i < 3; ++i)
    {
        for (uint j = 0; j < 2; ++j)
        {
            m_boundaries(i, j)->finalize();
            m_boundaries(i, j)->setAsUninitialized();
        }
    }
}


const string Site::info(int xr, int yr, int zr, string desc) const
{

    stringstream s_full;

    s_full << str();
    s_full << "[" << NX() << " x " << NY() << " x " << NZ() << "] * ";

    if (isActive())
    {
        s_full << "Active";

        const SoluteParticle *particle = associatedParticle();

        s_full << " " << particle->particleStateName();

        s_full << " * Neighbors: ";

        for (uint i = 0; i < Site::nNeighborsLimit(); ++i)
        {
            s_full << particle->nNeighbors(i) << " ";

        }


    }

    else
    {
        s_full << "Deactive";
    }


    s_full << "\n";

    uint _min = ParticleStates::nStates;

    ucube nN;
    nN.copy_size(m_levelMatrix);

    Site * currentSite;
    for (uint i = 0; i < m_neighborhoodLength; ++i)
    {
        for (uint j = 0; j < m_neighborhoodLength; ++j)
        {
            for (uint k = 0; k < m_neighborhoodLength; ++k)
            {

                currentSite = neighborhood(i, j, k);

                //BLOCKED
                if (currentSite == NULL)
                {
                    nN(i, j, k) = _min + 1;
                }

                //THIS SITE
                else if (currentSite == this)
                {
                    nN(i, j, k) = _min + 2;
                }

                //MARKED SITE
                else if ((i == Site::nNeighborsLimit() + xr) && (j == Site::nNeighborsLimit() + yr) && (k == Site::nNeighborsLimit() + zr))
                {
                    nN(i, j, k) = _min + 3;
                }

                //ACTIVE PARTICLE
                else if (currentSite->isActive())
                {
                    KMCDebugger_Assert(currentSite->associatedParticle(), !=, NULL);

                    nN(i, j, k) = currentSite->associatedParticle()->particleState();
                }

                else
                {
                    nN(i, j, k) = _min;
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
    numberSearchRepl(ParticleStates::surface,      ParticleStates::shortNames.at(ParticleStates::surface));
    numberSearchRepl(ParticleStates::solvant,     ParticleStates::shortNames.at(ParticleStates::solvant));

    numberSearchRepl(_min+0, ".");
    numberSearchRepl(_min+1, "#");
    numberSearchRepl(_min+2, isActive() ? associatedParticle()->particleStateShortName() + "^" : "D");
    numberSearchRepl(_min+3, desc);

    s_full << s;

    string full_string = s_full.str();

    return full_string;

}


void Site::setInitialNNeighborsLimit(const uint &nNeighborsLimit, bool check)
{

    if (nNeighborsLimit >= min(uvec({NX(), NY(), NZ()}))/2 && check)
    {
        cerr << "Neighbor reach must be lower than half the minimum box dimension to avoid sites directly affecting themselves." << endl;
        KMCSolver::exit();
    }


    m_nNeighborsLimit = nNeighborsLimit;

    m_neighborhoodLength = 2*m_nNeighborsLimit + 1;


    m_originTransformVector = linspace<ivec>(-(int)m_nNeighborsLimit, m_nNeighborsLimit, m_neighborhoodLength);


    m_levelMatrix.set_size(m_neighborhoodLength, m_neighborhoodLength, m_neighborhoodLength);

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

                m_levelMatrix(i, j, k) = getLevel(std::abs(m_originTransformVector(i)),
                                                  std::abs(m_originTransformVector(j)),
                                                  std::abs(m_originTransformVector(k)));
            }
        }
    }

    DiffusionReaction::setupPotential();

}



void Site::setInitialBoundaries(const umat &boundaryMatrix)
{

    KMCDebugger_AssertBool(m_boundaries.empty(), "Boundaries need to be clearedd before they can be initialized.");


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

            case Boundary::ConcentrationWall:
                m_boundaries(XYZ, orientation) = new ConcentrationWall(XYZ, orientation);

                break;

            default:

                cerr << "Unknown boundary type " << m_boundaryTypes(XYZ, orientation) << endl;
                KMCSolver::exit();

                break;
            }

        }

        if (!Boundary::isCompatible(m_boundaryTypes(XYZ, 0), m_boundaryTypes(XYZ, 1)))
        {
            cerr << "Mismatch in boundaries for " << XYZ << "'th dimension: " << m_boundaryTypes.t();
            KMCSolver::exit();
        }
    }

}

void Site::setInitialBoundaries(const int boundaryType)
{
    setInitialBoundaries(Boundary::allBoundariesAs(boundaryType));
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

const uint &Site::N(const uint i)
{
    return m_solver->N(i);
}




KMCSolver* Site::m_solver;

uint       Site::m_nNeighborsLimit = KMCSolver::UNSET_UINT;

uint       Site::m_neighborhoodLength = KMCSolver::UNSET_UINT;


ucube      Site::m_levelMatrix;

ivec       Site::m_originTransformVector;



field<Boundary*> Site::m_boundaries;

field<const Setting*> Site::m_boundaryConfigs;

umat Site::m_boundaryTypes;

uint Site::refCounter = 0;




ostream & operator << (ostream& os, const Site& ss)
{
    if (&ss == NULL)
    {
        os << "NULL";
    }

    else
    {
        os << ss.str();
    }

    return os;
}

