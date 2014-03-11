#pragma once

#include <armadillo>
#include <libconfig_utils/libconfig_utils.h>

using namespace arma;
using namespace libconfig;

namespace kMC
{

class KMCSolver;

class Boundary
{
public:

    Boundary(const uint dimension, const uint orientation);

    virtual ~Boundary();


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


    static bool isBlocked(const uint xi)
    {
        return xi == BLOCKED_COORDINATE;
    }


    static bool isCompatible(const int type1, const int type2, bool reverse = true);

    static void setupLocations(const uint x, const uint y, const uint z, uvec3 &loc);

    static void setMainSolver(KMCSolver * solver);

    static void resetAll()

    {
        solver = NULL;
    }



    enum BoundaryTypes
    {
        Periodic,
        Edge,
        Wall,
        ConsentrationWall
    };


private:

    static uint BLOCKED_COORDINATE;

    static KMCSolver * solver;

    const uint m_dimension;
    const uint m_orientation;


protected:

    static KMCSolver * mainSolver()
    {
        return solver;
    }

    const static uint & NX();

    const static uint & NY();

    const static uint & NZ();

    static uint N(const uint i);

    uint span() const
    {
        return N(m_orientation);
    }

    const uint & orientation() const
    {
        return m_orientation;
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
