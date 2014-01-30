#ifndef DIFFUSIONREACTION_H
#define DIFFUSIONREACTION_H

#include "../reaction.h"

#include <libconfig_utils/libconfig_utils.h>

class DiffusionReaction : public Reaction
{
public:


    DiffusionReaction(Site *destination);

    static void loadPotential(const Setting & setting);

    static const cube & potential() {
        return m_potential;
    }

private:

    double EspN = 2;
    double EspNN = 1;
    double beta = 1.0;
    double mu = 1;

    static cube m_potential;

    Site* destination;

    const uint & x1 () {
        return destination->x();
    }


    const uint & y1 () {
        return destination->y();
    }

    const uint & z1 () {
        return destination->z();
    }


    // Reaction interface
public:
    void calcRate();
    bool isActive();
    void execute();
};

#endif // DIFFUSIONREACTION_H
