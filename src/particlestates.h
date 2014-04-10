#pragma once

#include <vector>
#include <string>

using namespace std;

namespace kMC
{

struct ParticleStates
{
    enum AllStates
    {
        solvant,
        surface,
        crystal
    };

    const static vector<string> names;
    const static vector<string> shortNames;

    const static uint nStates;

};

}
