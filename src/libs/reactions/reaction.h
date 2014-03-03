#pragma once

#include <sys/types.h>
#include <sstream>
#include <set>

#include <libconfig_utils/libconfig_utils.h>

class Site;
class KMCSolver;

class Reaction
{
public:

    Reaction(string name = "Reaction");

    virtual ~Reaction();

    const string name;


    const static double UNSET_RATE;


    static void setSolverPtr(KMCSolver * solver);

    static void loadConfig(const Setting & setting);


    virtual void setDirectUpdateFlags(const Site * changedSite) = 0;

    virtual void setImplicitUpdateFlags()
    {
        m_updateFlags.push_back(defaultUpdateFlag);
    }

    void selectTriumphingUpdateFlag();

    virtual bool isNotBlocked() const = 0;

    virtual bool allowedAtSite() = 0;

    virtual void calcRate() = 0;

    virtual void execute() = 0;


    virtual const string info(int xr = 0, int yr = 0, int zr = 0, string desc = "X")  const;


    static void resetAll()
    {
        IDcount = 0;
    }

    const static double & linearRateScale()
    {
        return m_linearRateScale;
    }


    void setSite(Site* site)
    {
        m_reactionSite = site;
    }

    const vector<int> & updateFlags() const
    {
        return m_updateFlags;
    }

    const int & updateFlag() const
    {
        return m_updateFlag;
    }

    const uint & ID() const
    {
        return m_ID;
    }

    const double &  rate() const
    {
        return m_rate;
    }

    const uint & x() const;

    const uint & y() const;

    const uint & z() const;

    const Site * reactionSite() const
    {
        return m_reactionSite;
    }

    bool isType(const string name) const
    {
        return name.compare(this->name) == 0;
    }


    virtual string getFinalizingDebugMessage() const;

    virtual string getInfoSnippet() const
    {
        return "-";
    }

    bool operator == (const Reaction & other)
    {
        return this == &other;
    }

    const string str() const
    {
        stringstream s;
        s << name << "@(" << x() << ", " << y() << ", " << z() << ") [" << getInfoSnippet() << "]";
        return s.str();
    }


    //! Update flags are given in the order such that the minimum of the flag set is the
    //! triumphant flag.
    enum AllUpdateFlags
    {
        UNSET_UPDATE_FLAG = -1,
        defaultUpdateFlag = 0
    };

    friend class testBed;

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

    vector<int> m_updateFlags;
    int      m_updateFlag;


};

ostream & operator << (ostream& os, const Reaction& ss);

