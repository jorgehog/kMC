#pragma once

#include "../kmcsolver.h"

#include "../../ignis/include/ignis.h"

namespace kMC
{

class KMCParticles : public ignis::PositionHandler<uint>
{
public:

    KMCParticles(KMCSolver *solver) :
        m_solver(solver)
    {

    }

    KMCSolver *solver() const
    {
        return m_solver;
    }

private:

    KMCSolver *m_solver;



    // PositionHandler interface
public:
    uint count() const
    {
        return Site::nCrystals() + Site::nSolutionParticles();
    }

    uint operator ()(const uint n, const uint d) const
    {
        cout << "fixme " << n << " " << d << endl;
        return 0;
    }

    uint  &operator ()(const uint n, const uint d)
    {
        cout << "fixme " << n << " " << d << endl;
        uint * a = new uint(10);
        return *a;
    }
};

}
