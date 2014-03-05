#pragma once


#include <libconfig_utils/libconfig_utils.h>

using namespace libconfig;

namespace kMC
{

class KMCSolver;

class Boundary
{
public:
    Boundary();


    static void setMainSolver(KMCSolver * solver);

    static bool isBlocked(const int xi)
    {
        return xi != BLOCKED_COORDINATE;
    }

    static bool isCompatible(const int type1, const int type2, bool reverse = true);


    virtual uint transformCoordinate(const int xi) = 0;

    virtual void loadConfig(const Setting& setting)
    {
        (void) setting;
    }

    virtual void update() {}



    enum BoundaryTypes
    {
        Periodic,
        Wall,
        ConsentrationWall
    };


private:

    static int BLOCKED_COORDINATE;

    static uint m_NX;
    static uint m_NY;
    static uint m_NZ;

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

};

}
