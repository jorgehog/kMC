#pragma once

#include <sys/types.h>
#include <limits>
#include <sstream>
#include <set>
#include <functional>

#include <libconfig_utils/libconfig_utils.h>


namespace kMC
{

class Site;
class SoluteParticle;
class KMCSolver;

class Reaction
{
public:

    Reaction(SoluteParticle * reactant);

    virtual ~Reaction();

    virtual string name() const
    {
        return "Reaction";
    }

    virtual string type() const
    {
        return name();
    }

    static void setMainSolver(KMCSolver * m_solver);

    static void loadConfig(const Setting & setting);

    static void setBeta(const double beta);

    static void setLinearRateScale(const double linearRateScale)
    {
        m_linearRateScale = linearRateScale;
    }

    virtual void setDirectUpdateFlags(const SoluteParticle *changedReactant, const uint level)
    {

        (void) changedReactant;
        (void) level;

        registerUpdateFlag(defaultUpdateFlag);

    }

    virtual bool isAllowed() const = 0;

    virtual void execute() = 0;

    virtual void reset();

    virtual void registerBetaChange(const double newBeta)
    {
        (void) newBeta;
    }

    virtual double calcRate() = 0;

    virtual const string info(int xr = 0, int yr = 0, int zr = 0, string desc = "X")  const;

    template<typename T>
    void registerUpdateFlag(T flag)
    {
        static_assert(std::is_scalar<T>::value, "invalid update flag.");

        int iflag = static_cast<int>(flag);

        if (iflag < m_updateFlag)
        {
            m_updateFlag = iflag;
        }
    }

    void setLastUsedEnergy();

    void setRate()
    {
        setRate(calcRate());
    }

    //Do not use unless you know what you are doing.
    void changeRate(const double rate)
    {
        setRate(rate);
    }

    void forceNewRate(const double rate)
    {
        m_rate = rate;
    }

    void disable()
    {
        setRate(Reaction::UNSET_RATE);
    }

    void setAddress(const uint address)
    {
        m_address = address;
    }

    const uint & address() const
    {
        return m_address;
    }


    const static uint & nReactions()
    {
        return refCount;
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

    template<typename T>
    void forceUpdateFlag(const T flag)
    {
        static_assert(std::is_scalar<T>::value, "invalid update flag.");
        m_updateFlag = static_cast<int>(flag);
    }

    const double &  rate() const
    {
        return m_rate;
    }

    static const double & beta()
    {
        return m_beta;
    }

    const uint & x() const;

    const uint & y() const;

    const uint & z() const;

    bool hasVacantStatus() const;

    bool isType(const string name) const
    {
        return name.compare(this->type()) == 0;
    }

    bool hasName(const string name) const
    {
        return name.compare(this->name()) == 0;
    }

    SoluteParticle * reactant() const
    {
        return m_reactant;
    }

    virtual string getInfoSnippet() const
    {
        return "-";
    }

    string propertyString() const;

    bool operator == (const Reaction & other) const
    {
        return this == &other;
    }

    const string str() const
    {
        stringstream s;
        s << name() << "@(" << setw(3) << x() << "," << setw(3) << y() << "," << setw(3) << z() << ") [" << getInfoSnippet() << "]";
        return s.str();
    }


    //! Update flags are given in the order such that the minimum of the flag set is the
    //! triumphant flag.
    enum AllUpdateFlags
    {
        UNSET_UPDATE_FLAG = std::numeric_limits<int>::max(),
        defaultUpdateFlag = 0
    };

    static const double UNSET_RATE;

    static const double UNSET_ENERGY;

    static const uint UNSET_ADDRESS;

    const static uint & NX();

    const static uint & NY();

    const static uint & NZ();

    const static uint & _refCount()
    {
        return refCount;
    }

private:

    static KMCSolver* m_solver;

    static double m_beta;
    static double m_linearRateScale;

    static uint refCount;

    SoluteParticle* m_reactant = NULL;

    //Can this be in the reaction site?
    double m_lastUsedEnergy;

    double m_rate;

    int m_updateFlag;

    uint m_address;

protected:

    void setRate(const double rate);

    static KMCSolver * solver()
    {
        return m_solver;
    }

    template<typename T1, typename T2>
    string unsetIf(const T1 val, const T2 unsetVal) const
    {
        stringstream s;

        s << setw(5);

        if (val == unsetVal)
        {
            s << "UNSET";
        }

        else
        {
            s << setprecision(1) << fixed << val;
        }

        return s.str();

    }


};

}

ostream & operator << (ostream& os, const kMC::Reaction& ss);

