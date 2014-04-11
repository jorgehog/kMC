#pragma once

#include "particlestates.h"

#include <vector>
#include <set>
#include <sys/types.h>
#include <armadillo>
#include <assert.h>

#include <libconfig_utils/libconfig_utils.h>


using namespace arma;


namespace kMC
{


class KMCSolver;
class Reaction;
class DiffusionReaction;
class SoluteParticle;
class Boundary;

class Site
{
public:

    Site(uint _x, uint _y, uint _z);

    ~Site();

    /*
     * Static non-trivial functions
     */

    static void setMainSolver(KMCSolver* solver);

    static void loadConfig(const Setting & setting);


    static void initializeBoundaries();

    static void updateBoundaries();


    static uint shellSize(const uint level)
    {
        return std::pow(2*(level + 1) + 1, 3) - std::pow(2*level + 1, 3);
    }

    static constexpr uint closestShellSize = 26;

    static uint maxNeighbors();


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


    void introduceNeighborhood();


    bool hasNeighboring(const int state) const;

    uint countNeighboring(int state) const;


    SoluteParticle* associatedParticle() const
    {
        return m_associatedParticle;
    }


    void distanceTo(const Site * other, int &dx, int &dy, int &dz, bool absolutes = false) const;

    uint maxDistanceTo(const Site * other) const;


    void reset();

    void clearNeighborhood();


    const string info(int xr = 0, int yr = 0, int zr = 0, string desc = "X") const;


    void forEachNeighborDo(function<void (Site *)> applyFunction) const;

    void forEachNeighborDo_sendIndices(function<void (Site *, uint, uint, uint)> applyFunction) const;



    /*
     * Misc. trivial functions
     */

    static const KMCSolver * solver()
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


    static const Boundary * boundaries(const uint xyz, const uint loc)
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

    const uint & x() const
    {
        return m_x;
    }

    const uint & y() const
    {
        return m_y;
    }

    const uint & z() const
    {
        return m_z;
    }

    const uint & r(const uint i) const
    {
        switch (i) {
        case 0:
            return m_x;
            break;
        case 1:
            return m_y;
            break;
        case 2:
            return m_z;
            break;
        default:
            break;
        }

        return m_x;

    }



    Site* neighborhood(const uint x, const uint y, const uint z) const
    {
        return m_neighborhood[x][y][z];
    }

    bool operator == (const Site & other) const
    {
        return this == &other;
    }

    const string str() const
    {
        stringstream s;
        s << "Site@(" << x() << ", " << y() << ", " << z() << ")";
        return s.str();
    }

    const static uint & NX();

    const static uint & NY();

    const static uint & NZ();

private:

    static field<Boundary*> m_boundaries;

    static field<const Setting*> m_boundaryConfigs;

    static umat m_boundaryTypes;


    static uint m_nNeighborsLimit;

    static uint m_neighborhoodLength;


    static ucube m_levelMatrix;

    static ivec m_originTransformVector;

    static KMCSolver* m_solver;


    Site**** m_neighborhood;


    SoluteParticle *m_associatedParticle;


    const uint m_x;

    const uint m_y;

    const uint m_z;

    static uint refCounter;

};

}

ostream& operator<<(ostream& os, const kMC::Site& ss);
