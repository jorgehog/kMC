#pragma once


#include "../reaction.h"

#include <libconfig_utils/libconfig_utils.h>


class DiffusionReaction : public Reaction
{
public:


    DiffusionReaction(Site *destination);

    ~DiffusionReaction() {}

    double getSaddleEnergy();

    static void loadConfig(const Setting & setting);

    static void resetAll()
    {
        m_potential.reset();
    }


    static const cube & potential()
    {
        return m_potential;
    }

    const uint & xD () const
    {
        return destination->x();
    }

    const uint & yD () const
    {
        return destination->y();
    }

    const uint & zD () const
    {
        return destination->z();
    }


    //tmp
    double lastUsedE = 0;
    double lastUsedEsp = 0;
    static uint counter;
    static uint total;


    friend class testBed;


private:

    static double rPower;
    static double scale;

    static cube m_potential;

    Site* destination;


    // Reaction interface
public:
    void calcRate();

    bool isNotBlocked();

    void execute();

    void dumpInfo(int xr = 0, int yr = 0, int zr = 0);

    bool allowedAtSite();

    bool isAffectedByChangeIn(const Site *site) const;

    string getInfoSnippet() const
    {
        stringstream s;

        s << xD() << "," << yD() << "," << zD();

        return s.str();
    }

};
