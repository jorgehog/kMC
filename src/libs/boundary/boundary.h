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

    static bool isBlocked(const uint xi)
    {
        return xi != BLOCKED_COORDINATE;
    }



    virtual uint transformCoordinate(const uint xi) = 0;

    virtual void loadConfig(const Setting& setting)
    {
        (void) setting;
    }

    virtual void update() {}

private:

    static uint BLOCKED_COORDINATE;

    static uint m_NX;
    static uint m_NY;
    static uint m_NZ;

    static KMCSolver * m_mainSolver;

protected:

    static KMCSolver * mainSolver()
    {
        return m_mainSolver;
    }

};

}
