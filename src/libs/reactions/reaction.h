#pragma once

#include <sys/types.h>
#include <sstream>
#include <set>

#include <libconfig_utils/libconfig_utils.h>


namespace kMC
{


class Site;
class KMCSolver;

class Reaction
{
public:

    Reaction(Site * currentSite, const string name = "Reaction");

    virtual ~Reaction();

    const string name;


    void addUpdateFlag(int flag)
    {
        m_updateFlags.insert(flag);
    }

    void selectTriumphingUpdateFlag();

    static void setMainSolver(KMCSolver * m_solver);

    static void loadConfig(const Setting & setting);


    virtual void setDirectUpdateFlags(const Site * changedSite) = 0;

    virtual bool isAllowed() const = 0;

    virtual void calcRate() = 0;

    virtual void execute() = 0;


    virtual const string info(int xr = 0, int yr = 0, int zr = 0, string desc = "X")  const;


    static void clearAll()
    {
        m_IDCount = 0;
    }

    const uint & IDCount()
    {
        return m_IDCount;
    }

    const static double & linearRateScale()
    {
        return m_linearRateScale;
    }

    const double & lastUsedEnergy() const
    {
        return m_lastUsedEnergy;
    }

    const set<int> & updateFlags() const
    {
        return m_updateFlags;
    }

    const int & updateFlag() const
    {
        return m_updateFlag;
    }

    void forceUpdateFlag(int flag)
    {
        m_updateFlag = flag;
    }

    const uint & ID() const
    {
        return m_ID;
    }

    const double &  rate() const
    {
        return m_rate;
    }

    const double & beta() const
    {
        return m_beta;
    }

    const uint & x() const;

    const uint & y() const;

    const uint & z() const;

    bool isType(const string name) const
    {
        return name.compare(this->name) == 0;
    }

    const Site * getReactionSite() const
    {
        return m_reactionSite;
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

    static const double UNSET_RATE;
    static const double UNSET_ENERGY;

    const static uint & NX();

    const static uint & NY();

    const static uint & NZ();

private:

    static KMCSolver* m_solver;

    static double m_beta;
    static double m_linearRateScale;

    static uint m_IDCount;

    const uint m_ID;

    Site* m_reactionSite = NULL;

    double m_lastUsedEnergy;

    double m_rate;

    set<int> m_updateFlags;
    int      m_updateFlag;

protected:

    void setRate(const double rate);

    static KMCSolver * solver()
    {
        return m_solver;
    }

    Site * reactionSite() const
    {
        return m_reactionSite;
    }



};

}

ostream & operator << (ostream& os, const kMC::Reaction& ss);

