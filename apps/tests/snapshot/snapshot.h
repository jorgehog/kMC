#pragma once

#include <kMC>
#include <vector>

using namespace std;
using namespace kMC;

class SnapShot
{
public:

    SnapShot(KMCSolver* solver);

    static bool _switch;

    ucube siteBox;

    vector<double> allRates;
    vector<uvec> allreactions;

    uint timeWhenTaken;

    seed_type seed;

    bool operator==(const SnapShot& other) const;

};

ostream& operator<<(ostream& os, const SnapShot& ss);
