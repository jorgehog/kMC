#pragma once


#include "../site.h"
#include <sys/types.h>

#include <libconfig_utils/libconfig_utils.h>

class Site;
class KMCSolver;

class Reaction
{
public:

    Reaction(string name = "Reaction");

    virtual ~Reaction();

    const string name;

    static void setSolverPtr(KMCSolver * solver);

    static void loadConfig(const Setting & setting);

    virtual void setUpdateFlags(const Site * changedSite) = 0;

    virtual bool isNotBlocked() = 0;

    virtual bool allowedAtSite() = 0;

    virtual void calcRate() = 0;

    virtual void execute() = 0;


    virtual void dumpInfo(int xr = 0, int yr = 0, int zr = 0)  const;


    virtual void setSite(Site* site)
    {
        m_reactionSite = site;
    }

    static void resetAll()
    {
        IDcount = 0;
    }

    const uint & ID() const
    {
        return m_ID;
    }

    const double &  rate() const
    {
        return m_rate;
    }

    const static double & linearRateScale()
    {
        return m_linearRateScale;
    }

    const uint & x() const
    {
        return m_reactionSite->x();
    }

    const uint & y() const
    {
        return m_reactionSite->y();
    }

    const uint & z() const
    {
        return m_reactionSite->z();
    }

    const Site * reactionSite() const
    {
        return m_reactionSite;
    }

    const static double UNSET_RATE;

    virtual string getInfoSnippet() const
    {
        return "-";
    }

    bool operator == (const Reaction & other)
    {
        return this == &other;
    }


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

ostream & operator << (ostream& os, const Reaction& ss);

