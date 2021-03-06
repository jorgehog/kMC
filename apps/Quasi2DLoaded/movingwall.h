#pragma once

#include <kMC>
#include <set>


namespace kMC
{

class QuasiDiffusionSystem;
class QuasiDiffusionReaction;

class MovingWall : public KMCEvent
{
public:

    MovingWall(const double E0,
               const double sigma0,
               const double r0,
               const ivec &heighmap);

    const string numericDescription() const
    {
        stringstream s;
        s << "E0_" << m_E0
          << "_s0_" << m_s0
          << "_r0_" << m_r0;

        return s.str();

    }

    ~MovingWall();

    static double expSmallArg(double arg)
    {
        if (arg > 0.1 || arg < -0.1)
        {
            return exp(arg);
        }

        BADAssClose(arg, 0, 0.1, "Argument is not small.", [&arg] ()
        {
            BADAssSimpleDump(arg);
        });

        double arg2 = arg*arg;
        double arg4 = arg2*arg2;
        double approx = 1.0 + arg*(1 + 1.0/6*arg2 + 1.0/120*arg4) + 0.5*(arg2 + 1.0/12*arg4);

        BADAssClose(exp(arg), approx, 1E-5,
                    "Exponential approximation failed.", [&] ()
        {
            BADAssSimpleDump(arg, exp(arg), approx);
        });

        return approx;
    }

    void initialize();

    void execute()
    {

    }

    void reset();

    void markAsAffected(SoluteParticle *particle)
    {
        m_affectedParticles.insert(particle);
    }

    bool isAffected(SoluteParticle *particle)
    {
        return m_affectedParticles.find(particle) != m_affectedParticles.end();
    }

    const set<SoluteParticle *> affectedParticles() const
    {
        return m_affectedParticles;
    }

    const double &changeInHeight() const
    {
        return m_dh;
    }

    double localPressure(const uint site) const
    {
        if (m_onsetTime != 0) BADAssClose(m_localPressure(site), localPressureEvaluate(site), 1E-5);

        return m_localPressure(site);
    }

    double localPressureEvaluate(const uint site) const
    {
        return _pressureExpression(m_h - m_heighmap(site));
    }

    const double &height() const
    {
        return m_h;
    }

    const uint &length() const
    {
        return m_heighmap.n_elem;
    }

    uint span() const
    {
        return m_h - m_heighmap.min();
    }

    uint cavityVolume() const
    {
        uint V = 0;

        for (uint site = 0; site < length(); ++site)
        {
            V += floor(m_h - m_heighmap(site));
        }

        return V;
    }

    static double r0FromEs(const double h0, const double EsMax, const double EsInit)
    {
        return (h0 - 1.0)/std::log(EsMax/EsInit);
    }

    static double s0FromEs(const double h0, const double EsMax, const double EsInit)
    {
        using std::pow;

        return pow(pow(EsMax, h0)/EsInit, (1./(h0 - 1)));
    }

    double pressureEnergySum() const
    {
        double s = 0;

        for (uint site = 0; site < m_heighmap.size(); ++site)
        {
            BADAssClose(m_localPressure(site), localPressureEvaluate(site), 1E-4);

            s += m_localPressure(site);
        }

        return s;
    }

    void onHeightsSet()
    {
        recalculateAllPressures();
    }

    void _setSystem(const QuasiDiffusionSystem *system)
    {
        m_system = system;
    }

private:

    const QuasiDiffusionSystem *m_system;

    vector<vector<QuasiDiffusionReaction *>> m_pressureAffectedReactions;

    double m_h;
    double m_dh;
    double m_mPrev;

    uint m_changed;

    const double m_r0;
    const double m_s0;
    const double m_E0;

    const ivec &m_heighmap;
    vec m_localPressure;
    set<SoluteParticle *> m_affectedParticles;

    void _rescaleHeight();

    void _updatePressureRates();

    void _locateNewAffected();

    void _removeBlocked();

    void remakeUpdatedValues();

    void recalculateAllPressures();

    double _pressureExpression(const double heightDifference) const
    {
        return -m_s0*std::exp(-heightDifference/m_r0);
    }

};


}
