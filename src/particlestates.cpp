#include "particlestates.h"

using namespace kMC;


const vector<string> ParticleStates::names = {"crystal", "fixedcrystal", "solution", "surface"};

const vector<string> ParticleStates::shortNames = {"C", "F", "P", "S"};


int ParticleStates::equalAs(int state)
{
    switch (state) {
    case ParticleStates::crystal:
        return ParticleStates::fixedCrystal;
        break;
    default:
        return state;
        break;
    }
}
