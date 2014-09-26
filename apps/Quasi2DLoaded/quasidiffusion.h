#pragma once

#include <kMC>
#include "movingwall.h"

using namespace arma;

namespace kMC
{

class QuasiDiffusionReaction : public Reaction
{
public:

    QuasiDiffusionReaction(SoluteParticle *particle, ivec &heights, const double Eb, MovingWall &wallEvent) :
        Reaction(particle),
        m_heights(heights),
        m_Eb(Eb),
        m_rightSite(rightSite(1)),
        m_leftSite(leftSite(1)),
        m_wallEvent(wallEvent)
    {
        registerUpdateFlag(UpdateFlags::CALCULATE);
    }

    virtual string name() const
    {
        return "QuasiDiffusionReaction";
    }

    virtual string type() const
    {
        return "QuasiDiffusionReaction";
    }

    const string numericDescription() const
    {
        stringstream s;
        s << "Eb_" << Eb()
          << "_beta_" << beta()
          << "_" << wallEvent().numericDescription()
          << "_concentration_" << concentration();

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

    double wallHeight() const
    {
        if (!m_wallEvent.hasStarted())
        {
            return std::numeric_limits<double>::max();
        }

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
        return localEnergy();
    }

    virtual double prefactor() const
    {
        return 1.0;
    }

    virtual double calcRate() override final
    {
        if (updateFlag() != (int)UpdateFlags::CALCULATE && rate() != UNSET_RATE)
        {
            return rate();
        }

        double Ea = activationEnergy();

        if (Ea == 0)
        {
            return prefactor();
        }
        else
        {
            return prefactor()*exp(-beta()*Ea);
        }

    }

    double calcRateBruteForce() const
    {
        return prefactor()*exp(-beta()*activationEnergy());
    }

    virtual bool pressureAffected() const
    {
        return true;
    }

    const MovingWall &wallEvent() const
    {
        return m_wallEvent;
    }

    const double &concentration() const
    {
        return solver()->targetConcentration();
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
        if (!m_wallEvent.hasStarted() || m_wallEvent.cycle() == 0)
        {
            return 0;
        }

        return m_wallEvent.localPressure(site());
    }

    virtual const string info(int xr = 0, int yr = 0, int zr = 0, string desc = "X") const
    {
        (void) xr;
        (void) yr;
        (void) zr;
        (void) desc;

        stringstream s;

        s <<"site " << site() << " rs ls "<< rightSite() << " " << leftSite() << " height " << myHeight() << " dh " << heightDifference(leftSite()) << " " << heightDifference(rightSite()) << " rate " << rate() << " pressure " << localPressure() << endl;
        s << "\n" << Reaction::info() << endl;

        return s.str();
    }

    enum class UpdateFlags
    {
        CALCULATE = -1
    };

protected:


    void registerHeightChange(const uint site, const int change)
    {
        m_heights(site) += change;
    }

    void markAsAffected(SoluteParticle *particle)
    {
        if (m_wallEvent.hasStarted())
        {
            m_wallEvent.markAsAffected(particle);
        }

        for (Reaction *reaction : particle->reactions())
        {
            if (reaction->isAllowed())
            {
                reaction->registerUpdateFlag(UpdateFlags::CALCULATE);
            }
        }
    }

    void markAsAffected()
    {
        markAsAffected(reactant());
        markAsAffected(solver()->getSite(leftSite(), 0, 0)->associatedParticle());
        markAsAffected(solver()->getSite(rightSite(), 0, 0)->associatedParticle());
    }

    MovingWall &wallEvent()
    {
        return m_wallEvent;
    }

    const int &heights(const uint i) const
    {
        return m_heights(i);
    }

private:

    ivec &m_heights;

    const double m_Eb;

    const uint m_rightSite;
    const uint m_leftSite;

    MovingWall &m_wallEvent;


};

class LeftHop : public QuasiDiffusionReaction
{
    using QuasiDiffusionReaction::QuasiDiffusionReaction;

public:

    virtual string name() const
    {
        return "LeftHop";
    }

    bool isAllowed() const
    {
        return heights(leftSite()) < myHeight();
    }

    void execute()
    {
        registerHeightChange(site(), -1);
        registerHeightChange(leftSite(), +1);

        markAsAffected();
        markAsAffected(solver()->getSite(leftSite(2), 0, 0)->associatedParticle());
    }


    const string info(int xr = 0, int yr = 0, int zr = 0, string desc = "X") const
    {
        return "Lefthop " + QuasiDiffusionReaction::info(xr, yr, zr, desc);
    }

};

class RightHop : public QuasiDiffusionReaction
{
    using QuasiDiffusionReaction::QuasiDiffusionReaction;

public:

    virtual string name() const
    {
        return "RightHop";
    }

    bool isAllowed() const
    {
        return heights(rightSite()) < myHeight();
    }

    void execute()
    {
        registerHeightChange(site(), -1);
        registerHeightChange(rightSite(), +1);

        markAsAffected();
        markAsAffected(solver()->getSite(rightSite(2), 0, 0)->associatedParticle());
    }

    const string info(int xr = 0, int yr = 0, int zr = 0, string desc = "X") const
    {
        return "Righthop " + QuasiDiffusionReaction::info(xr, yr, zr, desc);
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

    using QuasiDiffusionReaction::QuasiDiffusionReaction;

public:

    virtual string name() const
    {
        return "Deposition";
    }

    bool pressureAffected() const
    {
        return false;
    }

    bool isAllowed() const
    {
        return myHeight() < floor(wallHeight());
    }

    void execute()
    {
        registerHeightChange(site(), +1);

        markAsAffected();
    }

    double activationEnergy() const = 0;

    double prefactor() const = 0;

    const string info(int xr = 0, int yr = 0, int zr = 0, string desc = "X") const
    {
        return "Deposition " + QuasiDiffusionReaction::info(xr, yr, zr, desc);
    }

};

class DepositionSpesifiedChemicalPotential : public Deposition
{
public:
    DepositionSpesifiedChemicalPotential(SoluteParticle *particle,
                                         ivec &heights,
                                         const double Eb,
                                         MovingWall &wallEvent,
                                         const double chemicalPotentialDifference) :
        Deposition(particle, heights, Eb, wallEvent),
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
        return concentration();
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
        if (heightDifference(rightSite()) > 0)
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
        return nVacantNeighbors()*concentration();
    }

};


class Dissolution : public QuasiDiffusionReaction
{
    using QuasiDiffusionReaction::QuasiDiffusionReaction;

public:

    virtual string name() const
    {
        return "Dissolution";
    }

    bool isAllowed() const
    {
        return concentration() != 1;
    }

    double activationEnergy() const
    {
        return localEnergy() + localPressure();
    }

    void execute()
    {
        registerHeightChange(site(), -1);

        markAsAffected();
    }

    const string info(int xr = 0, int yr = 0, int zr = 0, string desc = "X") const
    {
        return "Dissolution " + QuasiDiffusionReaction::info(xr, yr, zr, desc);
    }
};

}
