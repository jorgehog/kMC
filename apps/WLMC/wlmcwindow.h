#pragma once

#include <sys/types.h>

namespace kMC
{

class WLMCwindow
{
public:
    WLMCwindow(const uint minBin, const uint maxBin);

private:

    const uint m_minBin;
    const uint m_maxBin;

};

}
