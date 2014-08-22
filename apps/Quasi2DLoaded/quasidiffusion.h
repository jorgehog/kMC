#pragma once

#include <kMC>

using namespace kMC;
using namespace arma;

class QuasiDiffusionReaction : public Reaction
{
public:

    QuasiDiffusionReaction(SoluteParticle *particle, ivec &heights) :
        Reaction(particle),
        m_heights(heights),
        m_rightSite(rightSite(1)),
        m_leftSite(leftSite(1))
    {

    }

    static const string potentialString()
    {
        stringstream s;
        s << DiffusionReaction::linearRateScale()
          << "_" << DiffusionReaction::beta()
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

    const uint &site() const
    {
        return reactant()->x();
    }

    const uint &nSites() const
    {
        return solver()->NX();
    }

    static void initialize()
    {
        m_leftRightRates.set_size(4);

        const double &l = DiffusionReaction::linearRateScale();
        const double &b = DiffusionReaction::beta();
        const double &s = DiffusionReaction::strength();

        m_leftRightRates(0) = 0;
        m_leftRightRates(1) = exp(-b*(1 + 1*s));
        m_leftRightRates(2) = exp(-b*(1 + 2*s));
        m_leftRightRates(3) = exp(-b*(1 + 3*s));

        m_depositionRate = l;
    }

protected:

    ivec &m_heights;

    void queueAffected()
    {
        reactant()->markAsAffected();

        solver()->getSite(leftSite(), 0, 0)->associatedParticle()->markAsAffected();
        solver()->getSite(rightSite(), 0, 0)->associatedParticle()->markAsAffected();
    }

    static const double &leftRightRate(const uint n)
    {
        return m_leftRightRates(n);
    }

    static const double &depositionRate()
    {
        return m_depositionRate;
    }

private:

    const uint m_rightSite;
    const uint m_leftSite;

    static vec m_leftRightRates;
    static double m_depositionRate;

};

class LeftHop : public QuasiDiffusionReaction
{
    using QuasiDiffusionReaction::QuasiDiffusionReaction;

    // Reaction interface
public:
    bool isAllowed() const
    {
        return m_heights(leftSite()) < myHeight();
    }

    void calcRate()
    {
//        setRate(leftRightRate(nNeighbors())/heightDifference(leftSite()));
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
        return m_heights(rightSite()) < myHeight();
    }

    void calcRate()
    {
//        setRate(leftRightRate(nNeighbors())/heightDifference(rightSite()));
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

class Deposition : public QuasiDiffusionReaction
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
        setRate(depositionRate());
    }

    void execute()
    {
        m_heights(site())++;

        queueAffected();
    }
};

class Dissolution : public QuasiDiffusionReaction
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
        //skipped
        disable();
    }

    void execute()
    {
        m_heights(site())--;

        //Skipped desorbtion for now.

        queueAffected();
    }
};


vec QuasiDiffusionReaction::m_leftRightRates;
double QuasiDiffusionReaction::m_depositionRate;
