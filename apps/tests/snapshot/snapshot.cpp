#include "snapshot.h"

#include <sys/time.h>

SnapShot::SnapShot(KMCSolver *solver)
{

    timeWhenTaken = time(NULL);
    seed = Seed::initialSeed;

    siteBox.set_size(solver->NX(), solver->NY(), solver->NZ());

    Site* currentSite;

    for (uint i = 0; i < siteBox.n_rows; ++i)
    {
        for (uint j = 0; j < siteBox.n_cols; ++j)
        {
            for (uint k = 0; k < siteBox.n_slices; ++k)
            {
                currentSite = solver->getSite(i, j, k);
                siteBox(i, j, k) = currentSite->isActive();

                if (currentSite->isActive())
                {
                    for (Reaction * r : currentSite->getAssociatedParticle()->reactions())
                    {
                        allRates.push_back(r->rate());
                        allreactions.push_back({r->x(), r->y(), r->z(),
                                                ((DiffusionReaction*)r)->xD(),
                                                ((DiffusionReaction*)r)->yD(),
                                                ((DiffusionReaction*)r)->zD()});
                    }
                }

            }
        }
    }

}

ostream & operator << (ostream& os, const SnapShot& ss)
{

    (void) ss;

    os << (SnapShot::_switch ? ":-)" : ":-(");

    SnapShot::_switch = !SnapShot::_switch;

    return os;
}

bool SnapShot::operator==(const SnapShot &other) const
{

    bool equal = true;

    if (seed != other.seed)
    {
        cout << "mismatch in snapshot seeds: " << seed << " != " << other.seed << endl;
        return false;
    }

    for (uint l = 0; l < allRates.size(); ++l)
    {
        double diff = abs(allRates.at(l) - other.allRates.at(l));
        if (diff > 1E-10)
        {
            cout << "mismatch in snapshot rates: " << allRates.at(l) << " != " << other.allRates.at(l) << " for l = " << l << " diff: "<< setprecision(16) << diff << " (";

            for (uint rc : allreactions.at(l))
            {
                cout << rc << " ";
            }

            cout << ")" << endl;

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
                    cout << "mismatch in snapshot sites: " << siteBox(i, j, k) << " != " << other.siteBox(i, j, k) << "  for i, j, k = " << i << " " << j << " " << k << endl;
                    return false;
                }

            }
        }
    }

    return equal;

}


bool SnapShot::_switch = true;
