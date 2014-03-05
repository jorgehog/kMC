#pragma once


#include <libconfig_utils/libconfig_utils.h>

using namespace libconfig;

namespace kMC
{

class Boundary
{
public:
    Boundary();

    static void loadConfig(const Setting& setting);

    static bool isBlocked(const uint xi)
    {
        return xi != BLOCKED_COORDINATE;
    }

    virtual uint transformCoordinate(const uint xi) = 0;


private:

    static uint m_NX;
    static uint m_NY;
    static uint m_NZ;

    static uint BLOCKED_COORDINATE;

};

}
