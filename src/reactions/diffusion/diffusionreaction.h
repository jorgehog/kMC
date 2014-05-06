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

    double getSaddleEnergy();

    double getSaddleEnergyContributionFrom(const SoluteParticle *particle);

    double getSaddleEnergyContributionFromNeighborAt(const uint i, const uint j, const uint k, const uint s1, const uint s2);

    static umat::fixed<3, 2> makeSaddleOverlapMatrix(const ivec &relCoor);

    static void loadConfig(const Setting & setting);

    static void setupPotential();

    static void clearAll()
    {
        m_potential.reset();
        m_saddlePotential.reset_objects();
        m_saddlePotential.reset();
        m_neighborSetIntersectionPoints.reset_objects();
        m_neighborSetIntersectionPoints.reset();
    }

    static const double & potential(const uint & x, const uint & y, const uint & z, const uint speciesA = 0, const uint speciesB = 0)
    {
        return m_potential(x, y, z)(speciesA, speciesB);
    }

    static const field<mat> & potentialBox()
    {
        return m_potential;
    }

    Site *destinationSite() const;

    const int &path(const int i) const
    {
        return m_path[i];
    }

    const double & lastUsedEsp() const
    {
        return m_lastUsedEsp;
    }

    uint xD() const;

    uint yD() const;

    uint zD() const;


    //static setters
    static void setBetaChangeScaleFactor(const double factor)
    {
        m_betaChangeScaleFactor = factor;
    }

    static void setPotentialParameters(const vector<double> &rPowers, const double scale, bool setup = true);


    string getFinalizingDebugMessage() const;

private:

    static mat m_rPowers;

    static double m_scale;

    static double m_betaChangeScaleFactor;

    static field<mat> m_potential;
    static field<field<mat>> m_saddlePotential;

    static field<umat::fixed<3, 2> > m_neighborSetIntersectionPoints;


    double m_lastUsedEsp;

    enum SpecificUpdateFlags
    {
        updateKeepSaddle = 2
    };

    uint m_saddleFieldIndices[3];

    int m_path[3];

    // Reaction interface
public:

    void setDirectUpdateFlags(const SoluteParticle* changedReactant, const uint level);

    void calcRate();

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
