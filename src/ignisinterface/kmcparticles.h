#pragma once

#include <ignis/include/ignis.h>

namespace kMC
{

class KMCSolver;

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
    uint count() const;

    uint operator ()(const uint n, const uint d) const;
    uint &operator ()(const uint n, const uint d);
};

}
