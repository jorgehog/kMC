#pragma once

#include <kMC>

using namespace kMC;
using namespace arma;

class QuasiDiffusionReaction : public Reaction
{
public:

    QuasiDiffusionReaction(SoluteParticle *particle, ivec *heights) :
        Reaction(particle),
        m_heights(*heights),
        m_rightSite(rightSite(1)),
        m_leftSite(leftSite(1))
    {

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
            return 2;
        }

        else if (leftHug || rightHug)
        {
            return 1;
        }

        else
        {
            return 0;
        }

    }

    int myHeight() const
    {
        return m_heights(site());
    }

    const uint &site() const
    {
        return reactant()->x();
    }

    const uint &nSites() const
    {
        return solver()->NX();
    }

    static void initialize(const int maxHeight)
    {
        m_leftRightRates.set_size(3);

        m_leftRightRates(0) = DiffusionReaction::linearRateScale();
        m_leftRightRates(1) = DiffusionReaction::linearRateScale()*exp(-DiffusionReaction::beta()*DiffusionReaction::strength());
        m_leftRightRates(2) = DiffusionReaction::linearRateScale()*exp(-2*DiffusionReaction::beta()*DiffusionReaction::strength());

        m_maxHeight = maxHeight;
    }

protected:

    ivec &m_heights;

    void queueAffected()
    {
        reactant()->markAsAffected();

        solver()->getSite(leftSite(), 0, 0)->associatedParticle()->markAsAffected();
        solver()->getSite(rightSite(), 0, 0)->associatedParticle()->markAsAffected();
    }

    const double &leftRightRate(const uint n) const
    {
        return m_leftRightRates(n);
    }

    static int m_maxHeight;

private:

    const uint m_rightSite;
    const uint m_leftSite;

    static vec m_leftRightRates;

};

class LeftHop : public QuasiDiffusionReaction
{
    using QuasiDiffusionReaction::QuasiDiffusionReaction;

    // Reaction interface
public:
    bool isAllowed() const
    {
        return m_heights(leftSite()) == myHeight() - 1;
    }

    void calcRate()
    {
        setRate(leftRightRate(nNeighbors()));
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
        return m_heights(rightSite()) == myHeight() - 1;
    }

    void calcRate()
    {
        setRate(leftRightRate(nNeighbors()));
    }

    void execute()
    {
        m_heights(site())--;
        m_heights(rightSite())++;

        queueAffected();
        solver()->getSite(rightSite(2), 0, 0)->associatedParticle()->markAsAffected();
    }

};

class Adsorbtion : public QuasiDiffusionReaction
{
    using QuasiDiffusionReaction::QuasiDiffusionReaction;

    // Reaction interface
public:
    bool isAllowed() const
    {
        return true;
    }

    void calcRate()
    {
        setRate(1.0/leftRightRate(nNeighbors()));
    }

    void execute()
    {
        m_heights(site())++;

        if (myHeight() > m_maxHeight)
        {
            m_maxHeight = myHeight();
        }

        queueAffected();
    }
};

class Desorbtion : public QuasiDiffusionReaction
{
    using QuasiDiffusionReaction::QuasiDiffusionReaction;


    // Reaction interface
public:
    bool isAllowed() const
    {
        return true;
    }

    void calcRate()
    {
        setRate(leftRightRate(nNeighbors()));
    }

    void execute()
    {
        m_heights(site())--;

        //Skipped desorbtion for now.

        queueAffected();
    }
};

class LoadTracker : public KMCEvent
{
public:
};


vec QuasiDiffusionReaction::m_leftRightRates;
int QuasiDiffusionReaction::m_maxHeight;
