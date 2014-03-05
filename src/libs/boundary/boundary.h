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

private:

    static uint m_NX;
    static uint m_NY;
    static uint m_NZ;

};

}
