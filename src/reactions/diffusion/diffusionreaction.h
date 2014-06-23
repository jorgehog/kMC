#pragma once


#include "../reaction.h"

#include <armadillo>

#include <libconfig_utils/libconfig_utils.h>

using namespace arma;


namespace kMC
{


class DiffusionReaction : public Reaction
{
public:


    DiffusionReaction(SoluteParticle *reactant, int dx, int dy, int dz);

    ~DiffusionReaction();

    string name() const
    {
        return "DiffusionReaction";
    }

    static void loadConfig(const Setting & setting);

    static void setupPotential();

    static void clearAll()
    {
        m_potential.reset();
    }

    static const double & potential(const uint x, const uint y, const uint z, const uint speciesA = 0, const uint speciesB = 0)
    {
        return m_potential(x, y, z)(speciesA, speciesB);
    }

    static const field<mat> & potentialBox()
    {
        return m_potential;
    }

    static const double &rPower(const uint i = 0, const uint j = 0)
    {
        return m_rPowers(i, j);
    }

    static const double &strength(const uint i = 0, const uint j = 0)
    {
        return m_strengths(i, j);
    }

    Site *destinationSite() const;

    const int &path(const int i) const
    {
        return m_path[i];
    }

    const uint &pathIndex(const int i) const
    {
        return m_pathIndices[i];
    }

    const double &pathLength() const
    {
        return m_pathLengths(m_pathIndices[0], m_pathIndices[1], m_pathIndices[2]);
    }

    uint xD() const;

    uint yD() const;

    uint zD() const;


    //static setters
    static void setBetaChangeScaleFactor(const double factor)
    {
        m_betaChangeScaleFactor = factor;
    }

    static void setPotentialParameters(const vector<double> &rPowers,
                                       const vector<double> &strenghts);


    string getFinalizingDebugMessage() const;

private:

    static mat m_rPowers;

    static mat m_strengths;

    static double m_betaChangeScaleFactor;

    static field<mat> m_potential;

    static cube::fixed<3, 3, 3> m_pathLengths;

    int m_path[3];

    uint m_pathIndices[3];

    // Reaction interface
public:

    void execute();

    bool isAllowed() const;

    void reset();

    void registerBetaChange(const double newBeta);

    const string info(int xr = 0, int yr = 0, int zr = 0, string desc = "X") const;

    string getInfoSnippet() const
    {
        stringstream s;

        s << setw(2) << m_path[0] << "," << setw(2) << m_path[1] << "," << setw(2) << m_path[2];

        if (destinationSite() == NULL)
        {
            s << " (b)";
        }

        else
        {
            s << "    ";
        }

        return s.str();
    }

};

}
