#pragma once

#include "kmcsolver.h"
#include "particlestates.h"

#include <vector>
#include <sys/types.h>
#include <armadillo>

#include <BADAss/badass.h>
#include <libconfig_utils/libconfig_utils.h>


using namespace arma;


namespace kMC
{

class Reaction;
class DiffusionReaction;
class SoluteParticle;
class Boundary;

class Site
{
public:

    Site();

    ~Site();

    /*
     * Static non-trivial functions
     */

    static void setMainSolver(KMCSolver* solver);

    static void loadConfig(const Setting & setting);


    static void initializeBoundaries();

    static void updateBoundaries();

    static bool boundariesIsInitialized();


    static uint shellSize(const uint level)
    {
        return std::pow(2*(level + 1) + 1, 3) - std::pow(2*level + 1, 3);
    }

    static constexpr uint closestShellSize = 26;

    static uint maxNeighbors()
    {
        return std::pow(m_neighborhoodLength, 3) - 1;
    }


    /*
     *  Misc static property functions
     */

    static uint getLevel(uint i, uint j, uint k);

    static umat getCurrentCrystalBoxTopology();



    /*
     * Init / Reset / clear static implementations
     */

    static void setInitialNNeighborsLimit(const uint & nNeighborsLimit, bool check = true);

    static void setInitialBoundaries(const umat & boundaryMatrix);

    static void setInitialBoundaries(const int boundaryType);



    static void resetNNeighborsLimitTo(const uint & nNeighborsLimit, bool check = true);

    static void resetBoundariesTo(const umat & boundaryMatrix);

    static void resetBoundariesTo(const int boundaryType);



    static void clearAll();

    static void clearBoundaries();

    static void finalizeBoundaries();


    /*
     * Non-trivial functions
     */

    static Site *getNextNeighbor(const uint x, const uint y, const uint z, const int dr, const uint dim);

    static ivec3 getSurfaceNormal(const uint x, const uint y, const uint z, const Site *phantomSite = NULL);

    static int detectSurfaceOrientation(const uint x, const uint y, const uint z, const uint dim, const Site *phantomSite = NULL);

    static int checkDirectionForCrystals(const uint x, const uint y, const uint z, const uint dim, const int orientation, const Site *phantomSite = NULL);



    static bool hasNeighboring(const uint x, const uint y, const uint z, const int state);

    static uint countNeighboring(const uint x, const uint y, const uint z, const int state);


    SoluteParticle* associatedParticle() const
    {
        return m_associatedParticle;
    }


    static const string info(const uint x, const uint y, const uint z, int xr = 0, int yr = 0, int zr = 0, string desc = "X");


    static void distanceBetween(const uint x0, const uint y0, const uint z0, const uint x1, const uint y1, const uint z1, int &dx, int &dy, int &dz, bool absolutes = false);

    static uint maxDistanceBetween(const uint x0, const uint y0, const uint z0, const uint x1, const uint y1, const uint z1);

    static void forEachNeighborDo_sendPath(const uint x, const uint y, const uint z, function<void (Site *, int, int, int)> applyFunction);

    static void forEachNeighborDo(const uint x, const uint y, const uint z, function<void (Site *)> applyFunction)
    {
        forEachNeighborDo_sendPath(x, y, z, [&applyFunction] (Site *site, int dx, int dy, int dz)
        {
            (void)dx;
            (void)dy;
            (void)dz;
            applyFunction(site);
        });
    }

    static void forEachNeighborDo_sendIndices(const uint x, const uint y, const uint z, function<void (Site *, uint, uint, uint)> applyFunction);

    static void forShellDo(const uint x, const uint y, const uint z, const int shellNumber, function<void(Site *, int, int, int)> applyFunction);


    /*
     * Misc. trivial functions
     */

    static KMCSolver * solver()
    {
        return m_solver;
    }

    static const uint &boundaryTypes(const uint i, const uint j = 0)
    {
        return m_boundaryTypes(i, j);
    }

    static const uint &nNeighborsLimit()
    {
        return m_nNeighborsLimit;
    }

    static const uint &neighborhoodLength()
    {
        return m_neighborhoodLength;
    }

    static const uint &levelMatrix(const uint i, const uint j, const uint k)
    {
        return m_levelMatrix(i, j, k);
    }

    static uint originTransformVector(const uint i)
    {
        return m_originTransformVector(i);
    }

    static const ivec & originTransformVector()
    {
        return m_originTransformVector;
    }

    static const uvec & neighborhoodIndices()
    {
        return m_neighborhoodIndices;
    }

    static const Boundary *boundaries(const uint xyz, const uint loc)
    {
        return m_boundaries(xyz, loc);
    }

    static const field<Boundary*> & boundaryField()
    {
        return m_boundaries;
    }

    void associateWith(SoluteParticle *particle)
    {
        m_associatedParticle = particle;
    }

    void desociate()
    {
        m_associatedParticle = NULL;
    }

    bool isActive() const
    {
        return m_associatedParticle != NULL;
    }

    static Site *neighborhood(const int x, const int y, const int z, const int xr, const int yr, const int zr)
    {
        BADAss(x, >=, 0);
        BADAss(y, >=, 0);
        BADAss(z, >=, 0);

        BADAss(x, <, (int)NX());
        BADAss(y, <, (int)NY());
        BADAss(z, <, (int)NZ());

        return m_solver->getSite(xr + x, yr + y, zr + z);
    }

    static Site *neighborhood_fromIndex(const int x, const int y, const int z, const uint xr, const uint yr, const uint zr)
    {
        return neighborhood(x, y, z, m_originTransformVector(xr), m_originTransformVector(yr), m_originTransformVector(zr));
    }


    bool operator == (const Site & other) const
    {
        return this == &other;
    }

    static const string str(const uint x, const uint y, const uint z);

    static const uint & NX();

    static const uint & NY();

    static const uint & NZ();

    static const uint &N(const uint i);

    static const uint & _refCount()
    {
        return refCounter;
    }

private:

    static field<Boundary*> m_boundaries;

    static field<const Setting*> m_boundaryConfigs;

    static umat m_boundaryTypes;

    static uint m_nNeighborsLimit;

    static uint m_neighborhoodLength;


    static ucube m_levelMatrix;

    static ivec m_originTransformVector;

    static uvec m_neighborhoodIndices;


    static KMCSolver* m_solver;


    SoluteParticle *m_associatedParticle;


    static uint refCounter;

};

}

ostream& operator<<(ostream& os, const kMC::Site& ss);
