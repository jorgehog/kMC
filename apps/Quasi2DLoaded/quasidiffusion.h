#pragma once

#include <kMC>

#include "movingwall.h"
#include "quasidiffusionsystem.h"

using namespace arma;

namespace kMC
{

class QuasiDiffusionReaction : public Reaction
{
public:

    QuasiDiffusionReaction(SoluteParticle *particle, QuasiDiffusionSystem &system) :
        Reaction(particle),
        m_system(system),
        m_rightSite(rightSite(1)),
        m_leftSite(leftSite(1))
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
        s << "alpha_" << betaEb()
          << "mu_" << mu()
          << "_" << wallEvent().numericDescription()
          << "_concentration_" << concentration();

        return s.str();
    }

    const double &betaEb() const
    {
        return m_system.alpha();
    }

    const double &mu() const
    {
        return m_system.log2C();
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
        return m_system.leftSite(site(), n);
    }

    uint rightSite(const uint n) const
    {
        return m_system.rightSite(site(), n);
    }

    int myHeight() const
    {
        return heights(site());
    }

    int heightDifference(const uint other) const
    {
        return heights(site()) - heights(other);
    }

    double wallHeight() const
    {
        if (!wallEvent().hasStarted())
        {
            return std::numeric_limits<double>::max();
        }

        return wallEvent().height();
    }

    const uint &site() const
    {
        return x();
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
            return prefactor()*exp(-Ea);
        }

    }

    double calcRateBruteForce() const
    {
        return prefactor()*exp(-betaEb()*activationEnergy());
    }

    virtual bool pressureAffected() const
    {
        return true;
    }

    const MovingWall &wallEvent() const
    {
        return m_system.wallEvent();
    }

    const double &concentration() const
    {
        return solver()->targetConcentration();
    }

    double localEnergy(const uint n) const
    {
        return betaEb()*(n + localPressure() - 2) + m_system.log2C();
    }

    double localEnergy() const
    {
        return localEnergy(nNeighbors());
    }

    uint nNeighbors() const
    {
        return m_system.nNeighbors(leftSite(), rightSite(), site());
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
        m_system.registerHeightChange(site, change);
    }

    void markAsAffected(SoluteParticle *particle)
    {
        if (wallEvent().hasStarted())
        {
            m_system.markAsAffected(particle);
        }
        else
        {
            for (Reaction *reaction : particle->reactions())
            {
                if (reaction->isAllowed())
                {
                    reaction->registerUpdateFlag(UpdateFlags::CALCULATE);
                }
            }
        }
    }

    void markAsAffected()
    {
        markAsAffected(reactant());
        markAsAffected(solver()->getSite(leftSite(), 0, 0)->associatedParticle());
        markAsAffected(solver()->getSite(rightSite(), 0, 0)->associatedParticle());
    }

    const int &heights(const uint i) const
    {
        return m_system.heights(i);
    }

private:

    QuasiDiffusionSystem &m_system;

    const uint m_rightSite;
    const uint m_leftSite;

    double localPressure() const
    {
        if (!wallEvent().hasStarted())
        {
            return 0;
        }

        return wallEvent().localPressure(site());
    }

};

class LeftHop : public QuasiDiffusionReaction
{
    using QuasiDiffusionReaction::QuasiDiffusionReaction;

public:

    virtual string name() const
    {
        return "LeftHop";
    }

    virtual bool isAllowed() const
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

    virtual bool isAllowed() const
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

class LeftHopDownOnly : public LeftHop
{
    using LeftHop::LeftHop;

public:
    double activationEnergy() const override
    {
        return localEnergy();
    }

};

class RightHopDownOnly : public RightHop
{
    using RightHop::RightHop;

public:

    double activationEnergy() const override
    {
        return localEnergy();
    }

};

class RightHopIsotropic : public RightHop
{
    using RightHop::RightHop;

public:

    double activationEnergy() const override
    {
        return localEnergy();
    }


    // Reaction interface
public:
    bool isAllowed() const
    {
        return true;
    }
};

class LeftHopIsotropic : public LeftHop
{
    using LeftHop::LeftHop;

public:
    double activationEnergy() const override
    {
        return localEnergy();
    }

public:
    bool isAllowed() const
    {
        return true;
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

class DepositionMirrorImageArhenius : public Deposition
{
    using Deposition::Deposition;

public:

    static double promoteFactor(const double n)
    {
        return (2 - n/2);
    }

    double activationEnergy() const
    {
        return 0.0;
    }

    double prefactor() const
    {
        return promoteFactor(nNeighbors());
    }

};

class DepositionMirrorImageArheniusNoShadowing : public Deposition
{
    using Deposition::Deposition;

public:

    static double promoteFactor()
    {
        return 1.0;
    }

    double activationEnergy() const
    {
        return 0.0;
    }

    double prefactor() const
    {
        return promoteFactor();
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
        return true;
    }

    double activationEnergy() const
    {
        return localEnergy();
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
