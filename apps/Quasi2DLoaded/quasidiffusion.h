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

    const string numericDescription() const
    {
        stringstream s;
        s << "Eb_" << Eb()
          << "_beta_" << beta()
          << "_" << wallEvent().numericDescription();

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

    virtual double activationEnergy() const
    {
        return localEnergy(nNeighbors());
    }

    virtual double prefactor() const
    {
        return 1.0;
    }

    void calcRate()
    {
        double Ea = activationEnergy();

        if (Ea == 0)
        {
            setRate(prefactor());
        }
        else
        {
            setRate(prefactor()*exp(-beta()*activationEnergy()));
        }

    }


protected:

    ivec &m_heights;

    void queueAffected()
  {
      reactant()->markAsAffected();

      solver()->getSite(leftSite(), 0, 0)->associatedParticle()->markAsAffected();
      solver()->getSite(rightSite(), 0, 0)->associatedParticle()->markAsAffected();
  }

    const MovingWall &wallEvent() const
    {
        return m_wallEvent;
    }

    const double &Eb() const
    {
        return m_Eb;
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

    void execute()
    {
        m_heights(site())--;
        m_heights(leftSite())++;

        queueAffected();
        solver()->getSite(leftSite(2), 0, 0)->associatedParticle()->markAsAffected();
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

    void execute()
    {
        m_heights(site())--;
        m_heights(rightSite())++;

        queueAffected();
        solver()->getSite(rightSite(2), 0, 0)->associatedParticle()->markAsAffected();
    }

};

class LeftHopPressurized : public LeftHop
{
    using LeftHop::LeftHop;

public:
    double activationEnergy() const override
    {
        return localEnergy() + localPressure();
    }

};

class RightHopPressurized : public RightHop
{
    using RightHop::RightHop;

public:
    double activationEnergy() const override
    {
        return localEnergy() + localPressure();
    }

};


class Deposition : public QuasiDiffusionReaction
{
public:

    Deposition(SoluteParticle *particle,
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
        return myHeight() < floor(wallHeight()) && m_concentrationController.concentration() != 0;
    }

    void execute()
    {
        m_concentrationController.onParticleRemoval(site());
        m_heights(site())++;

        queueAffected();
    }

    double activationEnergy() const = 0;

    double prefactor() const = 0;


protected:

    const ConcentrationControl &concentrationController() const
    {
        return m_concentrationController;
    }

private:

    ConcentrationControl &m_concentrationController;

};

class DepositionSpesifiedChemicalPotential : public Deposition
{
public:
    DepositionSpesifiedChemicalPotential(SoluteParticle *particle,
                                         ivec &heights,
                                         const double Eb,
                                         const MovingWall &wallEvent,
                                         const double chemicalPotentialDifference,
                                         ConcentrationControl &concentrationController) :
        Deposition(particle, heights, Eb, wallEvent, concentrationController),
        m_chemicalPotentialDifference(chemicalPotentialDifference)
    {

    }

    double activationEnergy() const
    {
        return 1.0 + 2*Eb();
    }
    double prefactor() const
    {
        return exp(beta()*m_chemicalPotentialDifference);
    }

private:

    const double m_chemicalPotentialDifference;

};

class DepositionPurelyFromConcentration : public Deposition
{
    using Deposition::Deposition;

public:

    double activationEnergy() const
    {
        return 1.0 + 2*Eb();
    }

    double prefactor() const
    {
        return concentrationController().concentration();
    }

};

class DepositionMirrorImageArhenius : public Deposition
{
    using Deposition::Deposition;

public:

    uint nVacantNeighbors() const
    {
        uint nvn = 1;

        if (heightDifference(leftSite()) > 0)
        {
            nvn++;
        }
        if (heightDifference(rightSite() > 0))
        {
            nvn++;
        }

        return nvn;

    }

    double activationEnergy() const
    {
        return 1.0;
    }

    double prefactor() const
    {
        return nVacantNeighbors()*concentrationController().concentration();
    }


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

    double activationEnergy() const
    {
        return localEnergy() + localPressure();
    }

    void execute()
    {
        m_concentrationController.onParticleAddition(site());
        m_heights(site())--;

        queueAffected();
    }

private:

    ConcentrationControl &m_concentrationController;
};
