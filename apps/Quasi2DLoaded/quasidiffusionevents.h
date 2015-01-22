#pragma once

#include <kMC>

namespace kMC
{

class QuasiDiffusionSystem;
class QuasiDiffusionReaction;
class Dissolution;

class NNeighbors : public KMCEvent
{
public:

    NNeighbors(const QuasiDiffusionSystem &system) :
        KMCEvent("nNeighbors", "", true),
        m_system(system)
    {

    }

    void initialize()
    {
        m_sum = 0;
    }

    void execute();

    const double &localValue() const
    {
        return m_localValue;
    }

private:

    const QuasiDiffusionSystem &m_system;

    double m_sum;
    double m_localValue;

};

class EqConc : public KMCEvent
{
public:

    EqConc(const bool shadowing) :
        KMCEvent("EqConc", "", true, true),
        m_shadowing(shadowing),
        m_neighbours(0),
        m_dissolutionRate(0),
        m_totalTime(0)
    {
        setDependency(solver()->solverEvent());
    }

    virtual ~EqConc()
    {
        m_dissolutionReactions.clear();
    }

    void initialize();

    void execute();

    void reset();

    void restart();

    double dMu() const
    {
        return log(m_dMu);
    }

private:

    const double m_shadowing;

    double m_dMu;

    double m_neighbours;
    double m_dissolutionRate;

    vector<Dissolution*> m_dissolutionReactions;

    double m_totalTime;

    void update();

    void resetCounters()
    {
        m_totalTime = 0;

        m_neighbours = 0;
        m_dissolutionRate = 0;
    }


};

class ConcEquilibriator : public KMCEvent
{
public:

    ConcEquilibriator(QuasiDiffusionSystem &system,
                      EqConc &eqConcEvent,
                      const uint N = 100,
                      const double gCrit = 1E-5,
                      const double treshold = 1E-5) :
        KMCEvent("ConcEquilibriator", "", true, true),
        m_system(system),
        m_eqConcEvent(eqConcEvent),
        m_N(N),
        m_gCrit(gCrit),
        m_treshold(treshold),
        m_doAverage(false),
        m_averageMu(0),
        m_averageMu2Sum(0),
        m_averageMuCount(0),
        m_logMuShiftValues(m_N)
    {
        setDependency(m_eqConcEvent);
    }

    void initialize();

    void execute();

    double flatness()
    {
        double g = 0;

        for (uint i = 1; i < m_N - 1; ++i)
        {
            g += fabs((m_logMuShiftValues[i + 1] - m_logMuShiftValues[i - 1])/2.0);
        }

        return g/(m_N - 2);
    }

    const double &averageMu() const
    {
        return m_averageMu;
    }

    const double &error() const
    {
        return m_error;
    }

private:

    QuasiDiffusionSystem &m_system;

    ivec m_initialHeights;

    EqConc &m_eqConcEvent;

    vector<QuasiDiffusionReaction*> m_affectedReactions;

    const uint m_N;

    const double m_gCrit;
    const double m_treshold;

    bool m_doAverage;
    double m_averageMu;
    double m_averageMu2Sum;
    double m_error;
    uint m_averageMuCount;

    vector<double> m_logMuShiftValues;

    uint m_counter;

    void initiateNextConcentrationLevel();

    vector<double> m_shifts;
    vector<double> m_values;
    vector<double> m_swaps;

    double m_prevShift;

    void finalizeAverages();

};

class Cumulant : public KMCEvent
{
public:
    Cumulant(const QuasiDiffusionSystem &system) :
        KMCEvent("cumulant", "", true, true),
        m_system(system)
    {

    }

    void initialize()
    {
        m_cumulant = 0;
    }

    void execute();

    double cumulant() const;

private:

    double m_cumulant;
    const QuasiDiffusionSystem &m_system;

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

        for (uint dx = 0; dx < m_M; ++dx)
        {
            for (uint site = 0; site < m_M; ++site)
            {
                m_acf(dx) += (m_heightmap(site) - meanHeight)*(m_heightmap(site + dx) - meanHeight);
            }
        }

        if (cycle()%solver()->mainLattice()->outputSpacing() == 0)
        {
            (m_acf/((cycle()+1)*m_M)).eval().save(m_filename);
        }
    }

private:

    const ivec &m_heightmap;
    vec m_acfLocal;

    const uint m_M;
    vec m_acf;

    const string m_filename;

};


class SurfaceSize : public KMCEvent
{
public:

    SurfaceSize(const QuasiDiffusionSystem &system) :
        KMCEvent("SurfaceSize", "l0", true, true),
        m_system(system)
    {

    }

    void initialize()
    {
        m_sum = 0;
        m_T0 = T() - dt();
    }

    void execute();

    const double &localValue() const
    {
        return m_localValue;
    }

private:

    const QuasiDiffusionSystem& m_system;

    double m_sum;
    double m_localValue;

    double m_T0;

};

class SurfaceSizeLocal : public KMCEvent
{
public:

    SurfaceSizeLocal() :
        KMCEvent("SurfaceSizeLocal", "l0", true, true)
    {

    }

    void execute()
    {
        setValue(dependency<SurfaceSize>("SurfaceSize")->localValue());
    }

};



}
