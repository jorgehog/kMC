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


    const static double UNSET_RATE;


    static void setSolverPtr(KMCSolver * solver);

    static void loadConfig(const Setting & setting);

    void initialize();

    void update();

    bool isActive()
    {
        return m_siteReactionArrayIndex != UNSET_ARRAY_INDEX;
    }

    virtual bool isAffectedByChangeIn(const Site* site) const = 0;

    virtual void setUpdateFlags(const Site * changedSite, uint level) = 0;

    void getTriumphingUpdateFlag();

    virtual bool isNotBlocked() = 0;

    virtual bool allowedAtSite() = 0;

    virtual void calcRate() = 0;

    virtual void execute() = 0;


    virtual const string info(int xr = 0, int yr = 0, int zr = 0, string desc = ".")  const;


<<<<<<< HEAD
=======
    static void resetAll()
    {
        IDcount = 0;
    }

    const static double & linearRateScale()
    {
        return m_linearRateScale;
    }


>>>>>>> experimental2
    void setSite(Site* site)
    {
        m_reactionSite = site;
    }

    void clearUpdateFlags()
    {
        m_updateFlags.clear();
        assert(m_updateFlags.empty());
    }

    void setSiteReactionArrayIndex(uint index)
    {
        m_siteReactionArrayIndex = index;
    }

    const uint & ID() const
    {
        return m_ID;
    }

    const uint & siteReactionArrayIndex()
    {
        return m_siteReactionArrayIndex;
    }

    const double &  rate() const
    {
        return m_rate;
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

<<<<<<< HEAD
    const static double UNSET_RATE;
    const static uint   UNSET_ARRAY_INDEX;
=======
    bool isType(const string name) const
    {
        return name.compare(this->name) == 0;
    }

    void setUpdateFlag(int flag)
    {
        m_updateFlag = flag;
    }

    virtual string getFinalizingDebugMessage() const;
>>>>>>> experimental2

    virtual string getInfoSnippet() const
    {
        return "-";
    }

    bool operator == (const Reaction & other)
    {
        return this == &other;
    }

<<<<<<< HEAD
    string str() const
    {

        stringstream s;
        s << name << "@(" << x() << "," << y() << "," << z() << ") [" << getInfoSnippet() << "]";

        return s.str();

    }

=======
    const string str() const
    {
        stringstream s;
        s << name << "@(" << x() << "," << y() << "," << z() << ") [" << getInfoSnippet() << "]";
        return s.str();
    }


    //! Update flags are given in the order such that the minimum of the flag set is the
    //! triumphant flag.
    enum AllUpdateFlags
    {
        defaultUpdateFlag = 0,
        noUpdate = 100,
        UNSET_UPDATE_FLAG = 1337
    };

>>>>>>> experimental2
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
    uint m_siteReactionArrayIndex;


    double m_rate;

<<<<<<< HEAD
    void enable();

    void disable();
=======
    set<int> m_updateFlags;
    int      m_updateFlag;

>>>>>>> experimental2

};

ostream & operator << (ostream& os, const Reaction& ss);

