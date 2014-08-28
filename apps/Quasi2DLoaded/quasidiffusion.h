#pragma once

#include <kMC>

using namespace kMC;
using namespace arma;

class MovingWall : public KMCEvent
{
public:

    MovingWall(const double h0,
               const double r0,
               const double s0,
               const ivec &heighmap) :
        KMCEvent("MovingWall", "l0", true, true),
        m_h0(h0),
        m_h(h0),
        m_r0(r0),
        m_s0(s0),
        m_heighmap(heighmap)
    {

    }

    void execute()
    {
        _rescaleHeight();

//        double E = 0;
//        for (uint site = 0; site < m_heighmap.size(); ++site)
//        {
//            E += localPressure(site);
//        }

//        setValue(E/(m_heighmap.size()*m_s0*exp(-m_h0/m_r0)));
    }

    double localPressure(const uint site) const
    {
        return -m_s0*std::exp(-(m_h - m_heighmap(site))/m_r0);
    }

    const double &height() const
    {
        return m_h;
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

        setValue(floor(m_h) - m2);
    }



};

class QuasiDiffusionReaction : public Reaction
{
public:

    QuasiDiffusionReaction(SoluteParticle *particle, ivec &heights, const MovingWall &wallEvent) :
        Reaction(particle),
        m_heights(heights),
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

    static void initialize()
    {
        m_leftRightRates.set_size(4);

        m_leftRightRates(0) = 0;
        m_leftRightRates(1) = exp(-beta()*localEnergy(1));
        m_leftRightRates(2) = exp(-beta()*localEnergy(2));
        m_leftRightRates(3) = exp(-beta()*localEnergy(3));

        m_depositionRate = DiffusionReaction::linearRateScale();
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

    static double localEnergy(const uint n)
    {
        return 1 + n*DiffusionReaction::strength();
    }

    double localEnergy() const
    {
        return localEnergy(nNeighbors());
    }

    static const double &depositionRate()
    {
        return m_depositionRate;
    }

    double localPressure() const
    {
        return m_wallEvent.localPressure(site());
    }

private:

    const uint m_rightSite;
    const uint m_leftSite;

    static vec m_leftRightRates;
    static double m_depositionRate;

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
    using QuasiDiffusionReaction::QuasiDiffusionReaction;

    // Reaction interface
public:
    bool isAllowed() const
    {
        return myHeight() < floor(wallHeight());
    }

    void calcRate()
    {
        setRate(depositionRate());
    }

    void execute()
    {
        m_heights(site())++;
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
        setRate(depositionRate()/2);
    }

    void execute()
    {
        m_heights(site())--;
    }
};


vec QuasiDiffusionReaction::m_leftRightRates;
double QuasiDiffusionReaction::m_depositionRate;
