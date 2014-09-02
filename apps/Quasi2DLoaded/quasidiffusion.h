#pragma once

#include <kMC>
#include "quasidiffusionevents.h"

using namespace kMC;
using namespace arma;

class QuasiDiffusionReaction : public Reaction
{
public:

    QuasiDiffusionReaction(SoluteParticle *particle, ivec &heights, const double Eb, const MovingWall &wallEvent) :
        Reaction(particle),
        m_heights(heights),
        m_Eb(Eb),
        m_rightSite(rightSite(1)),
        m_leftSite(leftSite(1)),
        m_wallEvent(wallEvent)
    {

    }

    static const string potentialString()
    {
        stringstream s;
        s << DiffusionReaction::linearRateScale()
          << "_" << beta()
          << "_" << DiffusionReaction::strength();

        return s.str();
    }


    const uint &leftSite() const
    {
        return m_leftSite;
    }

    const uint &rightSite() const
    {
        return m_rightSite;
    }

    uint leftSite(const uint n) const
    {
        return (site() + nSites() - n)%nSites();
    }

    uint rightSite(const uint n) const
    {
        return (site() + n)%nSites();
    }

    bool connectedLeft() const
    {
        return m_heights(leftSite()) >= myHeight();
    }

    bool connectedRight() const
    {
        return m_heights(rightSite()) >= myHeight();
    }

    uint nNeighbors() const
    {
        bool leftHug = connectedLeft();
        bool rightHug = connectedRight();

        if (leftHug && rightHug)
        {
            return 3;
        }

        else if (leftHug || rightHug)
        {
            return 2;
        }

        else
        {
            return 1;
        }

    }

    int myHeight() const
    {
        return m_heights(site());
    }

    int heightDifference(const uint other) const
    {
        return m_heights(site()) - m_heights(other);
    }

    const double &wallHeight() const
    {
        return m_wallEvent.height();
    }

    const uint &site() const
    {
        return reactant()->x();
    }

    const uint &nSites() const
    {
        return solver()->NX();
    }

protected:

    ivec &m_heights;


    const double &Eb() const
    {
        return m_Eb;
    }

    double leftRightRate(const uint n) const
    {
        return exp(-beta()*localEnergy(n));
    }

    double localEnergy(const uint n) const
    {
        return 1 + n*m_Eb;
    }

    double localEnergy() const
    {
        return localEnergy(nNeighbors());
    }

    double localPressure() const
    {
        return m_wallEvent.localPressure(site());
    }


private:

    const double m_Eb;

    const uint m_rightSite;
    const uint m_leftSite;

    const MovingWall &m_wallEvent;

};

class LeftHop : public QuasiDiffusionReaction
{
    using QuasiDiffusionReaction::QuasiDiffusionReaction;

    // Reaction interface
public:
    bool isAllowed() const
    {
        return m_heights(leftSite()) < myHeight() && m_heights(leftSite()) < floor(wallHeight());
    }

    virtual void calcRate()
    {
        setRate(leftRightRate(nNeighbors()));
    }

    void execute()
    {
        m_heights(site())--;
        m_heights(leftSite())++;
    }


};

class RightHop : public QuasiDiffusionReaction
{
    using QuasiDiffusionReaction::QuasiDiffusionReaction;

    // Reaction interface
public:
    bool isAllowed() const
    {
        return m_heights(rightSite()) < myHeight() && m_heights(rightSite()) < floor(wallHeight());
    }

    virtual void calcRate()
    {
        setRate(leftRightRate(nNeighbors()));
    }

    void execute()
    {
        m_heights(site())--;
        m_heights(rightSite())++;
    }

};


class LeftHopScaled : public LeftHop
{
    using LeftHop::LeftHop;

    // Reaction interface
public:

    void calcRate()
    {
        setRate(leftRightRate(nNeighbors())/heightDifference(leftSite()));
    }

};

class RightHopScaled : public RightHop
{
    using RightHop::RightHop;

    // Reaction interface
public:

    void calcRate()
    {
        setRate(leftRightRate(nNeighbors())/heightDifference(rightSite()));
    }

};


class LeftHopPressurized : public LeftHop
{
    using LeftHop::LeftHop;

    // Reaction interface
public:

    void calcRate()
    {
        setRate(exp(-beta()*(localEnergy() + localPressure())));
    }

};

class RightHopPressurized : public RightHop
{
    using RightHop::RightHop;

    // Reaction interface
public:

    void calcRate()
    {
        setRate(exp(-beta()*(localEnergy() + localPressure())));
    }

};


class Deposition : public QuasiDiffusionReaction
{
public:

    Deposition(SoluteParticle *particle,
               ivec &heights,
               const double Eb,
               const MovingWall &wallEvent,
               const double chemicalPotentialDifference,
               ConcentrationControl &concentrationController) :
        QuasiDiffusionReaction(particle,
                               heights,
                               Eb,
                               wallEvent),
        m_chemicalPotentialDifference(chemicalPotentialDifference),
        m_concentrationController(concentrationController)
    {

    }


    // Reaction interface
public:
    bool isAllowed() const
    {
        return myHeight() < floor(wallHeight()) && m_concentrationController.nSolvants() != 0;
    }

    void calcRate()
    {
        double mu = 1.0;
        double prefactor = mu*m_concentrationController.nSolvants()*exp(beta()*m_chemicalPotentialDifference);

        double activatioEnergy = 1.0 + 4*m_concentrationController.concentration()*Eb();

        setRate(prefactor*exp(-beta()*activatioEnergy));
    }

    void execute()
    {
        m_concentrationController.registerPopulationChange(-1);
        m_heights(site())++;
    }

private:

    const double m_chemicalPotentialDifference;

    ConcentrationControl &m_concentrationController;
};



class Dissolution : public QuasiDiffusionReaction
{
public:

    Dissolution(SoluteParticle *particle,
                ivec &heights,
                const double Eb,
                const MovingWall &wallEvent,
                ConcentrationControl &concentrationController) :
    QuasiDiffusionReaction(particle,
                           heights,
                           Eb,
                           wallEvent),
    m_concentrationController(concentrationController)
    {

    }

    // Reaction interface
public:
    bool isAllowed() const
    {
        return m_concentrationController.concentration() != 1;
    }

    void calcRate()
    {
        setRate(exp(-beta()*(localEnergy() + localPressure())));
    }

    void execute()
    {
        m_concentrationController.registerPopulationChange(+1);
        m_heights(site())--;
    }

private:

    ConcentrationControl &m_concentrationController;
};
