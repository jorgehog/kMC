#pragma once

#include <kMC>

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
        m_r0(r0FromEs(h0, EsMax, EsInit)),
        m_s0(s0FromEs(h0, EsMax, EsInit)),
        m_heighmap(heighmap)
    {

    }

    void execute()
    {
        _rescaleHeight();
    }

    double localPressure(const uint site) const
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

    uint freeVolume() const
    {
        uint V = 0;

        for (uint site = 0; site < length(); ++site)
        {
            V += floor(m_h - m_heighmap(site));
        }

        return V - length();
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

    const double m_h0;
    double m_h;

    const double m_r0;
    const double m_s0;

    const ivec &m_heighmap;


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

        m_h = m_r0*std::log(m) + m_h0;

        setValue((m_h - m2)/m_h0);
    }



};


class ConcentrationControl : public KMCEvent
{
public:

    ConcentrationControl(const double cBoundary,
                         const double diffusivity,
                         const double dx,
                         const MovingWall &movingWallEvent) :
        KMCEvent("ConcentrationControl", "", true, true),
        m_boundaryConcentration(cBoundary),
        m_diffusivity(diffusivity),
        m_movingWallEvent(movingWallEvent),
        m_dx(dx),
        m_dxSquared(dx*dx)
    {

    }

    virtual ~ConcentrationControl()
    {

    }

    virtual void onConfinementGeometricChange() = 0;

    virtual void onParticleAddition(const uint x) = 0;

    virtual void onParticleRemoval(const uint x) = 0;

    virtual double concentration() const = 0;

    uint nSolvants() const
    {
        return concentration()*m_movingWallEvent.freeVolume();
    }

protected:

    void execute()
    {
        //diffusion etc. can be added here.

        diffuse();

        setValue(concentration());

    }

    const double m_boundaryConcentration;

    const double m_diffusivity;

    const MovingWall &m_movingWallEvent;

    const double m_dx;
    const double m_dxSquared;

private:


    virtual void diffuse() = 0;

};

class ConcentrationControl1D : public ConcentrationControl
{
public:

    ConcentrationControl1D(const double cBoundary,
                           const double diffusivity,
                           const uint nCells,
                           const double concentrationFieldLength,
                           const MovingWall &movingWallEvent) :
        ConcentrationControl(cBoundary, diffusivity, concentrationFieldLength/nCells, movingWallEvent),
        m_nCells(nCells),
        m_concentrations(m_nCells)
    {

    }

    virtual ~ConcentrationControl1D()
    {

    }

    void initialize()
    {
        m_concentrations.fill(m_boundaryConcentration);
    }

    void registerPopulationChange(int change)
    {

        double dC = 1.0/m_movingWallEvent.freeVolume();

        m_concentrations(0) += change*dC;

        BADAssBool(concentration() >= 0 && concentration() <= 1,
                   "Illegal concentration values initiated.", [&] ()
        {
            BADAssSimpleDump(change, concentration(), m_movingWallEvent.freeVolume());
        });

    }

    void onParticleAddition(const uint x)
    {
        (void) x;
        registerPopulationChange(+1);
    }
    void onParticleRemoval(const uint x)
    {
        (void) x;
        registerPopulationChange(-1);
    }

    double concentration() const
    {
        return m_concentrations(0);
    }

    void onConfinementGeometricChange()
    {

    }

private:

    const uint m_nCells;

    vec m_concentrations;



    void diffuse()
    {
        vec newConcentration(m_nCells);

        double fac = m_diffusivity*solver()->solverEvent()->lastTimeStep()/(m_dxSquared);


        newConcentration(m_nCells - 1) = m_concentrations(m_nCells - 1) + fac*(m_boundaryConcentration - 2*m_concentrations(m_nCells - 1) + m_concentrations(m_nCells - 2));

        for (int cell = m_nCells - 2; cell > 0; --cell)
        {
            newConcentration(cell) = m_concentrations(cell) + fac*(m_concentrations(cell + 1) - 2*m_concentrations(cell) + m_concentrations(cell - 1));
        }

        double boundaryCondition = m_concentrations(0);
        newConcentration(0) = m_concentrations(0) + fac*(m_concentrations(1) - 2*m_concentrations(0) + boundaryCondition);

        m_concentrations = newConcentration;

        if (nTimesExecuted() % MainLattice::saveFileSpacing() == 0)
        {
            m_concentrations.save(solver()->filePath() + "conc0.arma");
        }
    }

};

class ConcentrationControl3D : public ConcentrationControl
{
public:

    ConcentrationControl3D(const double cBoundary,
                           const double diffusivity,
                           const uint nCells,
                           const double dH,
                           const MovingWall &movingWallEvent) :
        ConcentrationControl(cBoundary, diffusivity, dH/nCells, movingWallEvent),
        m_dH(dH),
        m_nCells(nCells),
        m_concentrationField(movingWallEvent.length(), 2*nCells + movingWallEvent.height(), r())
    {

    }

    uint r() const
    {
        return m_movingWallEvent.span()/2 + m_nCells;
    }

    // Event interface
public:
    void initialize()
    {
        m_concentrationField.fill(m_boundaryConcentration);
    }

    // ConcentrationControl interface
public:
    void onParticleAddition(const uint x)
    {
        (void) x;
    }

    void onParticleRemoval(const uint x)
    {
        (void) x;
    }

    double concentration() const
    {
        return 0;
    }

    void onConfinementGeometricChange()
    {

    }

private:

    void diffuse()
    {

        uint yCentered, z;

        uint s = m_movingWallEvent.span();
        uint H = 2*m_nCells + s;

        for (uint y = 0; y < H; ++y)
        {
            yCentered = y - m_movingWallEvent.height() - s/2;

            z = sqrt((H*H)/4 - yCentered*yCentered);

            for (uint x = 0; x < m_movingWallEvent.length(); ++x)
            {
                m_concentrationField(x, y, z) = sqrt(pow(((int)y - (int)yCentered), 2) + z*z);
            }
        }

        m_concentrationField.save(solver()->filePath() + "conc3D.arma", raw_ascii);
        exit(1);
    }

    const double m_dH;

    const uint m_nCells;

    cube m_concentrationField;



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
        KMCEvent(),
        m_heighmap(heighmap),
        m_filename(solver()->filePath() + "heighmap.arma")
    {

    }

protected:

    void execute()
    {
        if (nTimesExecuted()%MainLattice::outputSpacing() == 0)
        {
            m_heighmap.save(m_filename);
        }
    }

private:

    const ivec &m_heighmap;

    const string m_filename;

};

}
