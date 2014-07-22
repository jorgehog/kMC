#pragma once

#include "window.h"

namespace WLMC
{

struct WindowParams
{
    uint m_lowerLimit;
    uint m_upperLimit;

    Window::OverlapTypes m_overlapType;

    bool m_allowSubwindowing;

    WindowParams(uint lowerLimit, uint upperLimit, Window::OverlapTypes overlapType) :
        m_lowerLimit(lowerLimit),
        m_upperLimit(upperLimit),
        m_overlapType(overlapType),
        m_allowSubwindowing(true)
    {

    }

};

}
