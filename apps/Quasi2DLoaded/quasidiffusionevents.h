#pragma once

#include <kMC>
#include <set>

class QuasiDiffusionReaction;

namespace kMC
{

class MovingWall : public KMCEvent
{
public:

    MovingWall(const double h0,
               const double EsMax,
               const double EsInit,
               const ivec &heighmap);

    const string numericDescription() const
    {
        stringstream s;
        s << "EsMax_" << m_EsMax
          << "_EsInit_" << m_EsInit
          << "_h0_" << m_h0;

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

    void execute();

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

    const double &localPressure(const uint site) const
    {
        return m_localPressure(site);
    }

    double localPressureEvaluate(const uint site) const
    {
        return -m_s0*std::exp(-(m_h - m_heighmap(site))/m_r0);
    }

    const double &initialHeight() const
    {
        return m_h0;
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


private:

    vector<vector<QuasiDiffusionReaction *>> m_pressureAffectedReactions;

    const double m_h0;
    double m_h;
    double m_dh;
    double m_mPrev;

    const double m_EsMax;
    const double m_EsInit;

    const double m_r0;
    const double m_s0;

    const ivec &m_heighmap;
    vec m_localPressure;
    set<SoluteParticle *> m_affectedParticles;

    void _rescaleHeight();

    void _updatePressureRates();

    void remakeUpdatedValues();

};


class HeightRMS : public KMCEvent
{
public:

    HeightRMS(const ivec &heightmap) :
        KMCEvent("heightRMS", "l0", true, true),
        m_heightmap(heightmap),
        m_L(heightmap.size())
    {

    }


protected:

    void execute()
    {

        const double &meanHeight = dependency("height")->value();

        double RMS = 0;

        for (uint i = 0; i < m_L; ++i)
        {
            RMS += (m_heightmap(i) - meanHeight)*(m_heightmap(i) - meanHeight);
        }

        RMS /= m_L;

        RMS = std::sqrt(RMS);

        setValue(RMS);

    }

private:

    const ivec &m_heightmap;
    const uint m_L;

};

class DumpHeighmap : public KMCEvent
{
public:

    DumpHeighmap(const ivec &heighmap) :
        KMCEvent("height", "", true, true),
        m_heighmap(heighmap),
        m_filename(solver()->filePath() + "heighmap.arma")
    {

    }

protected:

    void execute()
    {
        setValue(mean(conv_to<vec>::from(m_heighmap)));

        if (cycle()%solver()->mainLattice()->outputSpacing() == 0)
        {
            m_heighmap.save(m_filename);
        }
    }

private:

    const ivec &m_heighmap;

    const string m_filename;

};

class AutoCorrHeight : public KMCEvent
{
public:

    AutoCorrHeight(const ivec &heightmap) :
        KMCEvent("AutoCorrHeight"),
        m_heightmap(heightmap),
        m_M(heightmap.n_elem/2),
        m_acf(m_M),
        m_filename(solver()->filePath() + "acf.arma")
    {

    }

    const vec acf() const
    {
        return m_acf/cycle();
    }

    void execute()
    {
        const double &meanHeight = dependency("height")->value();

        vec acf = zeros(m_M);

        for (uint dx = 0; dx < m_M; ++dx)
        {
            for (uint site = 0; site < m_M; ++site)
            {
                acf(dx) += (m_heightmap(site) - meanHeight)*(m_heightmap(site + dx) - meanHeight);
            }
        }

        acf /= double(m_M);

        m_acf += acf;

        if (cycle()%solver()->mainLattice()->outputSpacing() == 0)
        {
            m_acf.save(m_filename);
        }
    }

private:

    const ivec &m_heightmap;

    const uint m_M;
    vec m_acf;

    const string m_filename;

};

}
