
#pragma once

#include <armadillo>

using namespace std;
using namespace arma;

namespace kMC
{

class KMCSolver;
class Site;

class Boundary
{
public:

    Boundary(const uint dimension, const uint orientation, const uint type);

    virtual ~Boundary();

    const uint type;

    virtual uint transformCoordinate(const int xi) const;

    virtual int getDistanceBetween(int x1, int x2)
    {
        return x1 - x2;
    }

    virtual void update() = 0;

    virtual void initialize() = 0;

    virtual void finalize() {}

    void setupBoundarySites();

    void distanceFromSite(const Site * site, int & dxi, bool abs = false);

    static bool isBlocked(const uint xi)
    {
        return (xi == BLOCKED_COORDINATE);
    }

    static bool isBlocked(const uint xi, const uint yi, const uint zi)
    {
        return isBlocked(xi) || isBlocked(yi) || isBlocked(zi);
    }

    static bool isCompatible(const int type1, const int type2, bool reverse = true);

    static umat allBoundariesAs(const int boundaryType)
    {
        return zeros<umat>(3, 2) + boundaryType;
    }

    static void setMainSolver(KMCSolver * m_solver);

    static void clearAll()
    {
        m_currentBoundaries.clear();

        m_solver = NULL;

    }



    enum BoundaryTypes
    {
        Periodic,
        Edge,
        ConcentrationWall
    };


    const static uint & NX();

    const static uint & NY();

    const static uint & NZ();

    static uint N(const uint i);


    static void setupCurrentBoundary(const uint x, const uint dim, const uint shift = 0);

    static void setupCurrentBoundaries(const uint x, const uint y, const uint z, const uint shift = 0);


    static const Boundary* currentBoundaries(const uint i)
    {
        return m_currentBoundaries.at(i);
    }

    const uint & orientation() const
    {
        return m_orientation;
    }

private:

    static uint BLOCKED_COORDINATE;

    static KMCSolver * m_solver;

    const uint m_dimension;

    const uint m_orientation;

    vector<Site*> m_boundarySites;


protected:

    static KMCSolver * solver()
    {
        return m_solver;
    }

    static uint getLocation(const uint xi, const uint dim, const uint shift);

    static vector<const Boundary*> m_currentBoundaries;

    vector<Site*> & boundarySites()
    {
        return m_boundarySites;
    }

    uint span() const
    {
        return N(m_dimension);
    }

    const uint & dimension() const
    {
        return m_dimension;
    }

    enum Dimensions
    {
        X,
        Y,
        Z
    };

    enum Orientations
    {
        Near,
        Far
    };


};

}
