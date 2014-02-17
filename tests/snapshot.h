#ifndef SNAPSHOT_H
#define SNAPSHOT_H

#include <kMC>
#include <vector>

using namespace std;

class SnapShot
{
public:

    SnapShot(KMCSolver* solver);

    ucube siteBox;

    vector<double> allRates;
    vector<uvec> allreactions;

    uint timeWhenTaken;

    bool operator==(const SnapShot& other) const;

};

ostream& operator<<(ostream& os, const SnapShot& ss);

#endif // SNAPSHOT_H
