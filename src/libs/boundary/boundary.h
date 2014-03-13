#pragma once

#include <armadillo>
#include <libconfig_utils/libconfig_utils.h>

using namespace arma;
using namespace libconfig;

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

    virtual uint transformCoordinate(const int xi) const
    {
        return (((xi >= (int)span()) || (xi < 0)) ? BLOCKED_COORDINATE : xi);
    }

    virtual int getDistanceBetween(int x1, int x2)
    {
        return x1 - x2;
    }

    virtual void loadConfig(const Setting& setting)
    {
        (void) setting;
    }

    virtual void update() = 0;

    virtual void initialize() = 0;

    void setupBoundarySites();

    void distanceFromSite(const Site * site, int & dxi, bool abs = false);


    static bool isBlocked(const uint xi)
    {
        return xi == BLOCKED_COORDINATE;
    }


    static bool isCompatible(const int type1, const int type2, bool reverse = true);


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
        Surface,
        ConcentrationWall
    };


    const static uint & NX();

    const static uint & NY();

    const static uint & NZ();

    static uint N(const uint i);


    static void setupCurrentBoundaries(const uint x, const uint y, const uint z);

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

    static void setupLocations(const uint x, const uint y, const uint z, uvec3 &loc);

    static vector<const Boundary*> m_currentBoundaries;

    const vector<Site*> & boundarySites() const
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

    enum Orientations
    {
        X,
        Y,
        Z
    };



};

}
