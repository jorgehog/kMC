#ifndef DIFFUSIONREACTION_H
#define DIFFUSIONREACTION_H

#include "../reaction.h"

class DiffusionReaction : public Reaction
{
public:

    static cube weights;

    DiffusionReaction(Site *destination);

private:

    double EspN = 2;
    double EspNN = 1;
    double beta = 4.0;
    double mu = 1;


    Site* destination;

    uint x1 () {
        return destination->x();
    }


    uint y1 () {
        return destination->y();
    }

    uint z1 () {
        return destination->z();
    }


    // Reaction interface
public:
    void calcRate();
    bool isActive();
    void execute();
};

#endif // DIFFUSIONREACTION_H
