#ifndef DIFFUSIONREACTION_H
#define DIFFUSIONREACTION_H

#include "../reaction.h"

#include <libconfig_utils/libconfig_utils.h>

class DiffusionReaction : public Reaction
{
public:


    DiffusionReaction(Site *destination);

    static void resetAll() {

        m_potential.reset();

    }

    static void loadPotential(const Setting & setting);

    static const cube & potential()
    {
        return m_potential;
    }

    const uint & x1 () {
        return destination->x();
    }

    const uint & y1 () {
        return destination->y();
    }

    const uint & z1 () {
        return destination->z();
    }

    double lastUsedE = 0;
    double lastUsedEsp = 0;

private:

    static double rPower;
    static double scale;

    static cube m_potential;

    Site* destination;

    double getSaddleEnergy();


    // Reaction interface
public:
    void calcRate();

    bool isNotBlocked();

    void execute();

    void dumpInfo(int xr = 0, int yr = 0, int zr = 0);

    void setupSiteDependencies() {
        destination->linkSiteDependency(reactionSite);
    }

    bool allowedAtSite();

    friend class testBed;
};

#endif // DIFFUSIONREACTION_H
