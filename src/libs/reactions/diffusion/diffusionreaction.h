#pragma once


#include "../reaction.h"

#include <libconfig_utils/libconfig_utils.h>


class DiffusionReaction : public Reaction
{
public:


    DiffusionReaction(Site *m_destinationSite);

    ~DiffusionReaction() {
        lastSetup.clear();
    }

    double getSaddleEnergy();

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

    const uint & xD () const
    {
        return m_destinationSite->x();
    }

    const uint & yD () const
    {
        return m_destinationSite->y();
    }

    const uint & zD () const
    {
        return m_destinationSite->z();
    }

    string getFinalizingDebugMessage() const;

    //tmp
    double lastUsedEnergy;
    double lastUsedEsp;
    set<const Site*> lastSetup;

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

    Site* m_destinationSite;

    enum SpecificUpdateFlags
    {
        updateKeepSaddle = 2
    };


    // Reaction interface
public:

    void setUpdateFlags(const Site * changedSite, uint level);

    void calcRate();

    bool isNotBlocked();

    void execute();

    const string info(int xr = 0, int yr = 0, int zr = 0, string desc = ".") const;

    bool allowedAtSite();

    string getInfoSnippet() const
    {
        stringstream s;

        s << xD() << "," << yD() << "," << zD();

        return s.str();
    }

};
