#pragma once

#include "diffusionreaction.h"

namespace kMC
{

class ArrheniusDiffusion : public DiffusionReaction
{
public:
    ArrheniusDiffusion(SoluteParticle *reactant, int dx, int dy, int dz);

    // DiffusionReaction interface
public:
    double calcRate();

    static void clearAll()
    {
        DiffusionReaction::clearAll();
    }
};


}
