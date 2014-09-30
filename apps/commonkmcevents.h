#pragma once

#include "../include/kMC"


namespace kMC
{

class RateChecker : public KMCEvent
{

public:

    RateChecker() :
        KMCEvent("rateChecker")
    {
        BADAssBool(KMCSolver::instanceActive(), "Create a KMCSolver before creating the rateChecker.");

        setDependency(KMCSolver::instance()->solverEvent());
    }

    void execute() {}

    void reset()
    {

        for (SoluteParticle *particle : solver()->particles())
        {
            for (Reaction* reaction : particle->reactions())
            {
                if (reaction->isAllowed())
                {
                    BADAss(reaction->address(), !=, Reaction::UNSET_ADDRESS);
                }
            }
        }

        for (Reaction *reaction : solver()->allPossibleReactions())
        {
            BADAssBool(reaction->isAllowed());

            double rate = reaction->rate();
            reaction->forceNewRate(Reaction::UNSET_RATE);

            double r2 = reaction->calcRate();
            BADAssClose(rate, r2, 1E-3, "error in rate updating.", [&] ()
            {
                reaction->forceNewRate(rate);
                cout << reaction->reactant()->info() << reaction->info() << endl;
                cout << "cr " << rate << " " << r2 << endl;
                KMCDebugger_DumpFullTrace(solver()->solverEvent()->selectedReaction()->info());
            });

            reaction->forceNewRate(rate);
        }

    }

};


class Concentration : public KMCEvent
{

public:

    Concentration(bool toFile = true) :
        KMCEvent("Concentration", "", true, toFile)
    {

    }

    void execute()
    {
        setValue(solver()->targetConcentration());
    }

};

class Sphericity : public KMCEvent
{
public:
    Sphericity() : KMCEvent("Sphericity", "", true, true) {}

protected:

    void execute()
    {
        uint A = 0;
        solver()->forEachSiteDo([&A] (uint x, uint y, uint z, Site* _site)
        {
            Site *site = _site;

            if (!site->isActive())
            {
                return;
            }

            if (site->associatedParticle()->isCrystal())
            {
                return;
            }

            if (!site->hasNeighboring(x, y, z, ParticleStates::crystal))
            {
                return;
            }

            uint nC = site->countNeighboring(x, y, z, ParticleStates::crystal);

            if (nC == 1)
            {
                A += 3;
            }

            else if (nC == 2 || nC == 3)
            {
                A += 2;
            }

            else
            {
                A += 1;
            }

        });

        setValue(pi3root*pow(6.0*(SoluteParticle::nCrystals() + A), 2./3)/A);
    }

private:

    static const double pi3root;

};

const double Sphericity::pi3root = pow(datum::pi, 1./3);

class tempChange : public KMCEvent
{
public:

    tempChange(const double T1, uint therm = 1) :
        KMCEvent("tempChange"),
        therm(therm),
        T1(T1)
    {

    }

    void initialize()
    {
        T0 = DiffusionReaction::beta();
        dT = (T1 - T0)/((eventLength()/(double)therm - 1));

    }

    void execute()
    {
        if (m_cycle%therm == 0)
        {
            DiffusionReaction::setBeta(T0 + dT*(m_cycle/therm));
        }
    }



private:

    uint therm;

    double T0;
    const double T1;

    double dT;


};


class TotalTime : public KMCEvent
{
public:

    TotalTime() : KMCEvent("Time", "s*", true, true) {}

protected:

    void execute()
    {
        setValue(solver()->solverEvent()->totalTime());
    }

};

class MeasureTemp : public KMCEvent
{
public:

    MeasureTemp() : KMCEvent("Temperature", "T*", true, true) {}

protected:

    void execute()
    {
        setValue(DiffusionReaction::beta());
    }

};

class AverageNeighbors : public KMCEvent
{
public:

    AverageNeighbors() : KMCEvent("avgN", "", true, true) {}

protected:

    void execute()
    {
        uint cN = 0;
        uint c  = 0;

        for (SoluteParticle *particle : solver()->particles())
        {
            cN += particle->nNeighbors();
            c++;
        }

        setValue(cN/(double(c)));
    }

};

class TotalEnergy : public KMCEvent
{
public:

    TotalEnergy() : KMCEvent("TotalEnergy", "E*", true, true) {}

protected:

    void execute()
    {
        setValue(SoluteParticle::totalEnergy());
    }

};


class Debug : public KMCEvent
{
public:

    Debug() : KMCEvent("debug", "", true) {}

protected:

    void execute()
    {
        setValue(registeredHandler(0, 0));
    }
};

class Sort : public KMCEvent
{
public:

    Sort() : KMCEvent() {}

protected:

    void execute()
    {
        if ((m_cycle+1)%250000 == 0)
        {
            solver()->sortReactionsByRate();
        }
    }

};

class CPUTime : public KMCEvent
{
public:

    CPUTime() : KMCEvent("CPUTime", "s", true, true) {}

    void initialize()
    {
        clock.tic();
    }

protected:


    void execute()
    {
        setValue(clock.toc());
    }

    arma::wall_clock clock;


};

class RandomInsertion : public KMCEvent
{
public:

    RandomInsertion() : KMCEvent() {}

    void initialize()
    {
        m_nPrev = SoluteParticle::nSolutionParticles();
    }

protected:

    void execute()
    {

        if (m_cycle%10 == 0)
        {
            if (SoluteParticle::nSolutionParticles() == 0)
            {
                solver()->insertRandomParticle();
            }
        }

        m_nPrev = SoluteParticle::nSolutionParticles();
    }

private:

    uint m_nPrev;

};


}
