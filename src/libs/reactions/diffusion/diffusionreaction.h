#pragma once


#include "../reaction.h"

#include <armadillo>

#include <libconfig_utils/libconfig_utils.h>

using namespace arma;

class DiffusionReaction : public Reaction
{
public:


    DiffusionReaction(Site *currentSite, Site *destinationSite);

    ~DiffusionReaction() {

    }

    static const double UNSET_ENERGY;

    double getSaddleEnergy();

    static umat::fixed<3, 2> getSaddleOverlapMatrix(const ivec &relCoor);

    static void loadConfig(const Setting & setting);

    static void resetAll()
    {
        counterAllRate = 0;
        counterEqSP = 0;
        totalSP = 0;
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

    const uint & xD () const;

    const uint & yD () const;

    const uint & zD () const;

    string getFinalizingDebugMessage() const;

    //tmp
    double lastUsedEnergy;
    double lastUsedEsp;

    static uint counterEqSP;
    static uint totalSP;

    static uint counterAllRate;

    static wall_clock timer;
    static double totalTime;


    friend class testBed;


private:

    static double rPower;
    static double scale;

    static cube m_potential;
    static field<cube> m_saddlePotential;


    Site* m_destinationSite;

    enum SpecificUpdateFlags
    {
        updateKeepSaddle = 2
    };

    umat::fixed<3, 2> neighborSetIntersectionPoints;
    ivec::fixed<3> path;
    vec::fixed<3> rSaddle;



    // Reaction interface
public:

    void setDirectUpdateFlags(const Site * changedSite);

    void calcRate();

    bool isNotBlocked() const;

    void execute();

    const string info(int xr = 0, int yr = 0, int zr = 0, string desc = "X") const;

    bool allowedAtSite();

    string getInfoSnippet() const
    {
        stringstream s;

        s << xD() << ", " << yD() << ", " << zD();

        return s.str();
    }

};
