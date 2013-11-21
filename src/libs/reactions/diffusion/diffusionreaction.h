#ifndef DIFFUSIONREACTION_H
#define DIFFUSIONREACTION_H

#include "../reaction.h"

class DiffusionReaction : public Reaction
{
public:
    DiffusionReaction(Site *destination);

private:

    double EspN = 1;
    double EspNN = 0;
    double temperature = 10;
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
