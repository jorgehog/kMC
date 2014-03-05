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

    Boundary();

    virtual ~Boundary();


    virtual uint transformCoordinate(const int xi) const = 0;

    virtual int getDistanceBetween(int x1, int x2)
    {
        return x1 - x2;
    }

    virtual void loadConfig(const Setting& setting)
    {
        (void) setting;
    }

    virtual void update() {}



    static bool isBlocked(const uint xi)
    {
        return xi == BLOCKED_COORDINATE;
    }


    static bool isCompatible(const int type1, const int type2, bool reverse = true);

    static void setupLocations(const uint x, const uint y, const uint z, uvec3 &loc);

    static void setMainSolver(KMCSolver * solver);

    static void reset()

    {
        m_NXYZ.clear();
        m_mainSolver = NULL;
    }



    enum BoundaryTypes
    {
        Periodic,
        Wall,
        ConsentrationWall
    };


private:

    static uint BLOCKED_COORDINATE;

    static uint m_NX;
    static uint m_NY;
    static uint m_NZ;

    static uvec m_NXYZ;

    static KMCSolver * m_mainSolver;

protected:

    static KMCSolver * mainSolver()
    {
        return m_mainSolver;
    }

    const uint & NX() const
    {
        return m_NX;
    }

    const uint & NY() const
    {
        return m_NY;
    }

    const uint & NZ() const
    {
        return m_NZ;
    }

    const uvec & NXYZ() const
    {
        return m_NXYZ;
    }

};

}
