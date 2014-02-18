#pragma once


#include "../reaction.h"

#include <libconfig_utils/libconfig_utils.h>


class DiffusionReaction : public Reaction
{
public:


    DiffusionReaction(Site *destination);

    ~DiffusionReaction() {}

    double getSaddleEnergy();

    static void loadConfig(const Setting & setting);

    static void resetAll()
    {
        m_potential.reset();
    }


    static const cube & potential()
    {
        return m_potential;
    }

    virtual void setSite(Site* site)
    {
        Reaction::setSite(site);

        int xT, yT, zT;

        site->distanceTo(destination, xT, yT, zT);

        site->addDiffusionReaction(this, 1 + xT, 1 + yT, 1 + zT);
    }

    const uint & xD ()
    {
        return destination->x();
    }

    const uint & yD ()
    {
        return destination->y();
    }

    const uint & zD ()
    {
        return destination->z();
    }

    //tmp
    double lastUsedE = 0;
    double lastUsedEsp = 0;
    static uint counter;
    static uint total;


    friend class testBed;


private:

    static double rPower;
    static double scale;

    static cube m_potential;

    Site* destination;


    // Reaction interface
public:
    void calcRate();

    bool isNotBlocked();

    void execute();

    void dumpInfo(int xr = 0, int yr = 0, int zr = 0);

    bool allowedAtSite();

};
