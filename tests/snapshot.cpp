#include "snapshot.h"

#include <sys/time.h>

SnapShot::SnapShot(KMCSolver *solver)
{

    timeWhenTaken = time(NULL);

    siteBox.set_size(solver->getNX(), solver->getNY(), solver->getNZ());

    Site* currentSite;

    for (uint i = 0; i < siteBox.n_rows; ++i)
    {
        for (uint j = 0; j < siteBox.n_cols; ++j)
        {
            for (uint k = 0; k < siteBox.n_slices; ++k)
            {
                currentSite = solver->getSites()[i][j][k];
                siteBox(i, j, k) = currentSite->isActive();

                for (Reaction * r : currentSite->siteReactions())
                {
                    allRates.push_back(r->rate());
                    allreactions.push_back({r->x(), r->y(), r->z()});
                }

            }
        }
    }

}

ostream & operator << (ostream& os, const SnapShot& ss)
{
    os << "snapshot@" << &ss;
    return os;
}

bool SnapShot::operator==(const SnapShot &other) const
{

    bool equal = true;

    for (uint l = 0; l < allRates.size(); ++l)
    {
        double diff = abs(allRates.at(l) - other.allRates.at(l));
        if (diff > 1E-15)
        {
            cout << "mismatch in snapshot rate calculation" << endl;
            cout << allRates.at(l) << " != " << other.allRates.at(l) << " for l = " << l << endl;
            cout << "diff: "<< setprecision(16) << diff << endl;

            cout << allreactions.at(l).t();

            equal = false;
            break;
        }
    }

    for (uint i = 0; i < siteBox.n_rows; ++i)
    {
        for (uint j = 0; j < siteBox.n_cols; ++j)
        {
            for (uint k = 0; k < siteBox.n_slices; ++k)
            {

                if (siteBox(i, j, k) != other.siteBox(i, j, k))
                {
                    cout << "mismatch in snapshot sites" << endl;
                    cout << siteBox(i, j, k) << " != " << other.siteBox(i, j, k);
                    cout << "  for i, j, k = " << i << " " << j << " " << k << endl;
                    return false;
                }

            }
        }
    }

    cout << "equal siteboxes" << endl;
    return equal;

}
