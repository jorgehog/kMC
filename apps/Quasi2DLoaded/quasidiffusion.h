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

//        queueAffected();
//        solver()->getSite(leftSite(2), 0, 0)->associatedParticle()->markAsAffected();

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

//        queueAffected();
//        solver()->getSite(rightSite(2), 0, 0)->associatedParticle()->markAsAffected();
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
        setRate(exp(-beta()*(localEnergy() + localPressure()))/heightDifference(leftSite()));
    }

};

class RightHopPressurized : public RightHop
{
    using RightHop::RightHop;

    // Reaction interface
public:

    void calcRate()
    {
        setRate(leftRightRate(exp(-beta()*(localEnergy() + localPressure())))/heightDifference(rightSite()));
    }

};


class Deposition : public QuasiDiffusionReaction
{
public:

    Deposition(SoluteParticle *particle,
               ivec &heights,
               const double Eb,
               const MovingWall &wallEvent,
               const double depositionRate) :
        QuasiDiffusionReaction(particle,
                               heights,
                               Eb,
                               wallEvent),
        m_depositionRate(depositionRate)
    {

    }


    // Reaction interface
public:
    bool isAllowed() const
    {
        return myHeight() < floor(wallHeight());
    }

    void calcRate()
    {
        setRate(m_depositionRate);
    }

    void execute()
    {
        m_heights(site())++;
    }

private:

    const double m_depositionRate;
};



class Dissolution : public QuasiDiffusionReaction
{
public:

    Dissolution(SoluteParticle *particle,
               ivec &heights,
               const double Eb,
               const MovingWall &wallEvent,
               const double dissolutionPrefactor) :
        QuasiDiffusionReaction(particle,
                               heights,
                               Eb,
                               wallEvent),
        m_dissolutionPrefactor(dissolutionPrefactor)
    {

    }

    // Reaction interface
public:
    bool isAllowed() const
    {
        return true;
    }

    void calcRate()
    {
        setRate(m_dissolutionPrefactor*(exp(-beta()*(localEnergy() + localPressure()))));
    }

    void execute()
    {
        m_heights(site())--;
    }

private:

    const double m_dissolutionPrefactor;
};
