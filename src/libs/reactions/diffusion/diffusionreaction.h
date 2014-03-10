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


    DiffusionReaction(Site *currentSite, Site *destinationSite);

    ~DiffusionReaction();


    double getSaddleEnergy();

    double getSaddleEnergyContributionFrom(const Site* site);

    double getSaddleEnergyContributionFromNeighborAt(const uint &i, const uint &j, const uint &k);


    static uint separation()
    {
        return m_separation;
    }

    static umat::fixed<3, 2> getSaddleOverlapMatrix(const ivec &relCoor);

    static void loadConfig(const Setting & setting);

    static void resetAll()
    {
        m_potential.reset();
    }


    static const double & potential(const uint & x, const uint & y, const uint & z)
    {
        return m_potential(x, y, z);
    }

    static const cube & potentialBox()
    {
        return m_potential;
    }

    const Site* destinationSite() const
    {
        return m_destinationSite;
    }

    const double & lastUsedEsp() const
    {
        return m_lastUsedEsp;
    }

    const uint & xD () const;

    const uint & yD () const;

    const uint & zD () const;

    string getFinalizingDebugMessage() const;


private:

    static double rPower;
    static double scale;

    static uint m_separation;

    static cube m_potential;
    static field<cube> m_saddlePotential;


    double m_lastUsedEsp;

    Site* m_destinationSite = NULL;

    enum SpecificUpdateFlags
    {
        updateKeepSaddle = 2
    };

    umat::fixed<3, 2> neighborSetIntersectionPoints;
    ivec::fixed<3> path;
    uvec::fixed<3> saddleFieldIndices;


    bool allowedGivenNotBlocked() const;

    // Reaction interface
public:

    void setDirectUpdateFlags(const Site * changedSite);

    void calcRate();

    void execute();

    bool isAllowed() const;

    const string info(int xr = 0, int yr = 0, int zr = 0, string desc = "X") const;

    string getInfoSnippet() const
    {
        stringstream s;

        s << xD() << ", " << yD() << ", " << zD();

        return s.str();
    }

};

}
