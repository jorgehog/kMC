#pragma once

#include <armadillo>

namespace kMC
{

class MovingWall;

class ConcentrationControl
{
public:

    ConcentrationControl(const double cBoundary,
                         const double diffusivity,
                         const double dx) :
        m_boundaryConcentration(cBoundary),
        m_diffusivity(diffusivity),
        m_dx(dx),
        m_dxSquared(dx*dx)
    {

    }


    virtual ~ConcentrationControl()
    {

    }

    virtual void initialize() = 0;

    virtual void onParticleAddition(const uint x) = 0;

    virtual void onParticleRemoval(const uint x) = 0;

    virtual double concentration() const = 0;

    uint nSolvants() const;

    virtual void diffuse(const double dt) = 0;

    void setMovingWallEvent(const MovingWall *movingWall)
    {
        m_movingWall = movingWall;
    }

protected:

    const double m_boundaryConcentration;

    const double m_diffusivity;

    const double m_dx;
    const double m_dxSquared;

    const MovingWall *m_movingWall;

};

class NoControl : public ConcentrationControl
{
public:
    NoControl(const double concentration) :
        ConcentrationControl(concentration, 1, 1)
    {

    }



    // ConcentrationControl interface
public:
    void initialize()
    {
    }
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
        return m_boundaryConcentration;
    }
    void diffuse(const double dt)
    {
        (void) dt;
    }
};


class ConcentrationControl1D : public ConcentrationControl
{
public:

    ConcentrationControl1D(const double cBoundary,
                           const double diffusivity,
                           const uint nCells,
                           const double concentrationFieldLength) :
        ConcentrationControl(cBoundary, diffusivity, concentrationFieldLength/nCells),
        m_nCells(nCells),
        m_concentrations(m_nCells)
    {

    }

    virtual ~ConcentrationControl1D()
    {

    }

    virtual void initialize() override final
    {
        m_concentrations.fill(m_boundaryConcentration);
    }

    void registerPopulationChange(int change);

    void onParticleAddition(const uint x) override final
    {
        (void) x;
        registerPopulationChange(+1);
    }
    void onParticleRemoval(const uint x) override final
    {
        (void) x;
        registerPopulationChange(-1);
    }

    double concentration() const override final
    {
        return m_concentrations(0);
    }

private:

    const uint m_nCells;

    arma::vec m_concentrations;



    virtual void diffuse(const double dt) override final;

};


class ConcentrationControl3D : public ConcentrationControl
{
public:

    ConcentrationControl3D(const double cBoundary,
                           const double diffusivity,
                           const uint nCells,
                           const double dH);

    virtual ~ConcentrationControl3D()
    {

    }

    uint r() const;

    virtual void initialize() override final;

    // ConcentrationControl interface
public:
    void onParticleAddition(const uint x) override final
    {
        (void) x;
    }

    void onParticleRemoval(const uint x) override final
    {
        (void) x;
    }

    double concentration() const override final
    {
        return 0;
    }

private:

    virtual void diffuse(const double dt) override final;

    const double m_dH;

    const uint m_nCells;

    arma::cube m_concentrationField;
};




}
