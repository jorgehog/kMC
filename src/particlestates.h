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
        crystal,
        fixedCrystal,
        solution,
        surface,
        any
    };

    const static vector<string> names;
    const static vector<string> shortNames;

    static int equalAs(int state);

};

}
