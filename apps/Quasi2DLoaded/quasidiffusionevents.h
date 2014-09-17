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
               const ivec &heighmap) :
        KMCEvent("MovingWall", "h0", true, true),
        m_h0(h0),
        m_h(h0),
        m_EsMax(EsMax),
        m_EsInit(EsInit),
        m_r0(r0FromEs(h0, EsMax, EsInit)),
        m_s0(s0FromEs(h0, EsMax, EsInit)),
        m_heighmap(heighmap),
        m_localPressure(heighmap.size())
    {

    }

    const string numericDescription() const
    {
        stringstream s;
        s << "EsMax_" << m_EsMax
          << "_EsInit_" << m_EsInit
          << "_h0_" << m_h0;

        return s.str();

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

    vector<QuasiDiffusionReaction *> m_pressureAffectedReactions;

    const double m_h0;
    double m_h;
    double m_dh;

    const double m_EsMax;
    const double m_EsInit;

    const double m_r0;
    const double m_s0;

    const ivec &m_heighmap;
    vec m_localPressure;
    set<SoluteParticle *> m_affectedParticles;

    void _rescaleHeight()
    {
        double m = 0;
        double m2 = 0;
        for (uint i = 0; i < m_heighmap.size(); ++i)
        {
            m += exp(m_heighmap(i)/m_r0);
            m2 += m_heighmap(i);
        }

        m /= m_heighmap.size();
        m2 /= m_heighmap.size();

        double hPrev = m_h;

        m_h = m_r0*std::log(m) + m_h0;

        m_dh = m_h - hPrev;

        setValue((m_h - m2)/m_h0);
    }

    void _updatePressureRates();

};


class heightRMS : public KMCEvent
{
public:

    heightRMS(const ivec &heightmap) :
        KMCEvent("heightRMS", "l0", true, true),
        m_heightmap(heightmap),
        m_L(heightmap.size())
    {

    }


protected:

    void execute()
    {

        double meanHeight = sum(m_heightmap)/double(m_heightmap.size());

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

}
