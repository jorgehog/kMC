#pragma once

#include <kMC>

namespace kMC
{

class Dissolution;

class EquilibriumConcentrationEstimator : public KMCEvent
{
public:

    EquilibriumConcentrationEstimator() :
        KMCEvent("EquilibriumConcentrationEstimator", "", true, true)
    {
        setDependency(solver()->solverEvent());
    }

    virtual ~EquilibriumConcentrationEstimator()
    {
        m_dissolutionReactions.clear();
    }

    void initialize();

    void execute();

    void reset();

    void restart()
    {
        m_eqConc = 0;

        m_neighbours = 0;
        m_dissolutionRate = 0;
    }

private:

    double m_eqConc;

    double m_neighbours;
    double m_dissolutionRate;

    double m_expFac;

    vector<const Dissolution*> m_dissolutionReactions;

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
