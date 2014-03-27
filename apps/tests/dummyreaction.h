#pragma once

#include <kMC>

namespace kMC
{


class DummyReaction : public Reaction
{
public:
    DummyReaction(const uint initAddress) :
        Reaction(NULL),
        initAddress(initAddress),
        allowed(true)
    {

    }

    // Reaction interface
public:
    void setDirectUpdateFlags(const Site *changedSite)
    {
        (void) changedSite;
    }

    bool isAllowed() const
    {
        return allowed;
    }

    bool isAllowedAndActive() const
    {
        return allowed;
    }

    void calcRate()
    {
        setRate(1.0);
    }

    void execute()
    {

    }

    uint initAddress;
    bool allowed;

};


}
