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

    double getSaddleEnergyContributionFromNeighborAt(const uint &i, const uint &j, const uint &k);

    static umat::fixed<3, 2> makeSaddleOverlapMatrix(const ivec &relCoor);

    static void loadConfig(const Setting & setting);

    static void setupPotential();

    static void clearAll()
    {
        m_potential.reset();
        m_saddlePotential.reset_objects();
        m_saddlePotential.reset();
        neighborSetIntersectionPoints.reset_objects();
        neighborSetIntersectionPoints.reset();
    }

    static const double & potential(const uint & x, const uint & y, const uint & z)
    {
        return m_potential(x, y, z);
    }

    static const cube & potentialBox()
    {
        return m_potential;
    }

    SoluteParticle *destination() const;

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

    static void setPotentialParameters(const double rPower, const double scale, bool setup = true)
    {

        m_rPower = rPower;

        m_scale = scale;

        if (setup)
        {
            setupPotential();
        }

    }


    string getFinalizingDebugMessage() const;

private:

    static double m_rPower;

    static double m_scale;

    static double m_betaChangeScaleFactor;

    static cube m_potential;
    static field<cube>  m_saddlePotential;
    static field<umat::fixed<3, 2> > neighborSetIntersectionPoints;


    double m_lastUsedEsp;

    enum SpecificUpdateFlags
    {
        updateKeepSaddle = 2
    };

    uint saddleFieldIndices[3];

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

        if (destination() == NULL)
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
