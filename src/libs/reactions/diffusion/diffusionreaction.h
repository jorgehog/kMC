#ifndef DIFFUSIONREACTION_H
#define DIFFUSIONREACTION_H

#include "../reaction.h"

#include <libconfig_utils/libconfig_utils.h>

class DiffusionReaction : public Reaction
{
public:


    DiffusionReaction(Site *destination);

    static void resetAll() {

        m_saddleTransformVector.reset();
        m_potential.reset();

    }

    static void loadPotential(const Setting & setting);

    static const cube & potential() {
        return m_potential;
    }

private:

    static double rPower;
    static double scale;

    static uint saddleCutoff;

    static cube m_potential;
    static ivec m_saddleTransformVector;

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

    double getSaddleEnergy();


    // Reaction interface
public:
    void calcRate();
    bool isActive();
    void execute();

    friend class testBed;
};

#endif // DIFFUSIONREACTION_H
