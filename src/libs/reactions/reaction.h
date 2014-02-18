#pragma once


#include "../site.h"
#include <sys/types.h>

#include <libconfig_utils/libconfig_utils.h>

class Site;
class KMCSolver;

class Reaction
{
public:

    Reaction();

    virtual ~Reaction();


    static void setSolverPtr(KMCSolver * solver);

    static void loadConfig(const Setting & setting);


    virtual bool isNotBlocked() = 0;

    virtual bool allowedAtSite()
    {
        return true;
    }

    virtual void calcRate() = 0;

    virtual void execute() = 0;


    virtual void dumpInfo(int xr = 0, int yr = 0, int zr = 0);


    virtual void setSite(Site* site)
    {
        m_reactionSite = site;
    }

    static void resetAll()
    {
        IDcount = 0;
    }

    const uint & ID()
    {
        return m_ID;
    }

    const double &  rate()
    {
        return m_rate;
    }

    const static double & linearRateScale()
    {
        return m_linearRateScale;
    }

    const uint & x()
    {
        return m_reactionSite->x();
    }

    const uint & y()
    {
        return m_reactionSite->y();
    }

    const uint & z()
    {
        return m_reactionSite->z();
    }

    const Site * reactionSite()
    {
        return m_reactionSite;
    }

    const static double UNSET_RATE;



protected:

    static KMCSolver* mainSolver;

    static uint NX;
    static uint NY;
    static uint NZ;

    static double beta;
    static double m_linearRateScale;

    static uint IDcount;


    uint m_ID;

    Site* m_reactionSite;

    double m_rate;

};
