#pragma once

#include <sys/types.h>
#include <climits>
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

    Reaction(Site * currentSite);

    virtual ~Reaction();

    static const string name;

    void registerUpdateFlag(int flag)
    {
        if (flag < m_updateFlag)
        {
            m_updateFlag = flag;
        }
    }

    void selectTriumphingUpdateFlag();

    static void setMainSolver(KMCSolver * m_solver);

    static void loadConfig(const Setting & setting);

    static void setBeta(const double beta)
    {
        m_beta = beta;
    }

    static void setLinearRateScale(const double linearRateScale)
    {
        m_linearRateScale = linearRateScale;
    }

    static void clearAll()
    {
        m_IDCount = 0;
    }

    virtual void setDirectUpdateFlags(const Site * changedSite) = 0;

    virtual bool isAllowed() const = 0;

    virtual void calcRate() = 0;

    virtual void execute() = 0;

    virtual void reset();


    virtual const string info(int xr = 0, int yr = 0, int zr = 0, string desc = "X")  const;



    const static uint & IDCount()
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

    const int & updateFlag() const
    {
        return m_updateFlag;
    }

    void resetUpdateFlag()
    {
        m_updateFlag = UNSET_UPDATE_FLAG;
    }

    void forceUpdateFlag(int flag)
    {
        m_updateFlag = flag;
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
        UNSET_UPDATE_FLAG = INT_MAX,
        defaultUpdateFlag = 0
    };

    static constexpr double UNSET_RATE = -1337;

    static constexpr double UNSET_ENERGY = -13371337;

    const static uint & NX();

    const static uint & NY();

    const static uint & NZ();

private:

    static KMCSolver* m_solver;

    static double m_beta;
    static double m_linearRateScale;

    static uint m_IDCount;

    Site* m_reactionSite = NULL;

    double m_lastUsedEnergy;

    double m_rate;

    int m_updateFlag;

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

