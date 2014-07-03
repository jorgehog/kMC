#pragma once


#include "site.h"

#include "reactions/reaction.h"

#include "debugger/debugger.h"

#include "RNG/kMCRNG.h"

#include "ignisinterface/kmcevent.h"

#include <sys/types.h>

#include <armadillo>

#include <limits>

#include <libconfig_utils/libconfig_utils.h>


using namespace arma;


class lammpswriter;

namespace kMC
{

const uint UNSET_UINT = std::numeric_limits<uint>::max();

class DumpFile;
class SolverEvent;

class KMCSolver
{
public:


    KMCSolver(const Setting & root);

    KMCSolver();

    ~KMCSolver();

    void reset();



    void mainloop()
    {
        m_mainLattice->eventLoop();
    }

    void addEvent(Event<uint> &event)
    {
        m_mainLattice->addEvent(event);
    }

    void addEvent(Event<uint> *event)
    {
        m_mainLattice->addEvent(event);
    }

    void addSubLattice(lattice &subLattice)
    {
        m_mainLattice->addSubField(subLattice);
    }

    void addSubLattice(lattice *subLattice)
    {
        m_mainLattice->addSubField(subLattice);
    }


    void setupMainLattice();


    bool spawnParticle(SoluteParticle *particle, const uint x, const uint y, const uint z, bool checkIfLegal);

    SoluteParticle *forceSpawnParticle(const uint x, const uint y, const uint z, const uint species = 0, const bool sticky = false);

    void despawnParticle(SoluteParticle *particle);

    const Reaction *executeRandomReaction();

    void initializeCrystal(const double relativeSeedSize, const uint species = 0, const bool sticky = false);

    void initializeSolutionBath(const uint species = 0, const bool sticky = false);

    void initializeLayers(const uint height, const uint start = 0, const uint species = 0, const bool sticky = false);

    void initializeFromXYZ(string path, uint frame);

    void initializeFromLAMMPS(string path, uint frame);

    void insertRandomParticle(const uint species = 0, const bool sticky = false, bool checkIfLegal = false);

    void forceRandomPosition(SoluteParticle *particle, bool checkIfLegal = false);

    void rotateSystem(const double yaw, const double pitch, const double roll);

    void sortParticles(function<bool(SoluteParticle*, SoluteParticle*)> comp)
    {
        std::sort(m_particles.begin(), m_particles.end(), comp);
    }


    void initializeParticles();


    void forEachSiteDo(function<void(uint, uint, uint, Site *)> applyFunction) const;


    void getRateVariables();

    uint getReactionChoice(const double R) const
    {
        return binarySearchForInterval(R, m_accuAllRates);
    }


    Site* getSite(const int i, const int j, const int k) const
    {
        return sites[i + m_boundaryPadding][j + m_boundaryPadding][k + m_boundaryPadding];
    }

    MainLattice *mainLattice() const
    {
        return m_mainLattice;
    }

    const vector<SoluteParticle*> & particles() const
    {
        return m_particles;
    }


    SoluteParticle *particle(const uint n) const
    {
        return m_particles.at(n);
    }

    const uint &NX () const
    {
        return m_NX;
    }

    const uint &NY () const
    {
        return m_NY;
    }

    const uint &NZ () const
    {
        return m_NZ;
    }

    const uint &N(const uint i) const
    {
        return m_N(i);
    }

    const uvec3 NVec() const
    {
        return m_N;
    }

    const uint &volume() const
    {
        return mainLattice()->volume;
    }

    const vector<double> & accuAllRates() const
    {
        return m_accuAllRates;
    }

    const vector<Reaction*> & allPossibleReactions() const
    {
        return m_allPossibleReactions;
    }

    const vector<uint> & availableReactionSlots() const
    {
        return m_availableReactionSlots;
    }

    const uint & cycle() const;

    const double & kTot() const
    {
        return m_kTot;
    }

    const double & targetConcentration()
    {
        return m_targetConcentration;
    }


    //Set functions

    void setBoxSize(const uvec3 boxSize, bool check = true)
    {
        setBoxSize(boxSize(0), boxSize(1), boxSize(2), check);
    }

    void setBoxSize(const uint NX, const uint NY, const uint NZ, bool check = true);


    void setNumberOfCycles(const uint nCycles)
    {
        MainLattice::nCycles = nCycles;
    }

    void setCyclesPerOutput(const uint cyclesPerOutput)
    {
        MainLattice::nCyclesPerOutput = cyclesPerOutput;
    }


    void setTargetConcentration(const double concentration)
    {
        if (concentration > 1 || concentration < 0)
        {
            exit("invalid concentration set,");
        }

        m_targetConcentration = concentration;
    }

    void setRNGSeed(uint seedState = Seed::fromTime, int defaultSeed = 0);


    static void exit(const string s = "Exit failure.")
    {
        throw std::runtime_error(s);
    }


    void clearParticles();


    void registerReactionChange(Reaction * reaction, const double &newRate);

    void reshuffleReactions();

    void swapReactionAddresses(const uint dest, const uint orig);

    void postReactionShuffleCleanup(const uint nVacancies);

    void updateAccuAllRateElements(const uint from, const uint to, const double value);

    void remakeAccuAllRates();

    double prevAccuAllRatesValue(const uint address) const
    {
        return address == 0 ? 0 : m_accuAllRates.at(address - 1);
    }

    bool isEmptyAddress(const uint address) const;

    bool isRegisteredParticle(SoluteParticle *particle) const;

    bool isPossibleReaction(Reaction *reaction) const;

    string getReactionVectorDebugMessage();

    void dumpXYZ(const uint n);

    void dumpLAMMPS(const uint n);


    static double minRateThreshold()
    {
        //        return max((*std::min_element(m_allPossibleReactions.begin(),
        //                                  m_allPossibleReactions.end(),
        //                                  [] (const Reaction *r1, const Reaction *r2) {return r1->rate() < r2->rate();}))->rate()/2,
        //                   1E-8);

        return Reaction::linearRateScale()*1E-8;

    }


    static void enableDumpXYZ(bool state)
    {
        m_dumpXYZ = state;
    }

    static void enableDumpLAMMPS(bool state)
    {
        m_dumpLAMMPS = state;
    }


    void initializeSites();

    void clearSites();

    const SolverEvent *solverEvent() const
    {
        return m_solverEvent;
    }

    void resetLastReaction();

    void sortReactionsByRate();

    static uint binarySearchForInterval(const double target, const vector<double> & intervals);


    void setFilepath(const string filepath)
    {
        m_filepath = filepath.empty() ? filepath : filepath + "/";
    }

    const string & filePath() const
    {
        return m_filepath;
    }

    enum class DiffusionTypes
    {
        None,
        TST,
        Arrhenius
    };

    const DiffusionTypes &diffusionType() const
    {
        return m_diffusionType;
    }

    void setDiffusionType(const DiffusionTypes type)
    {
        if (m_diffusionType == type)
        {
            return;
        }

        KMCDebugger_AssertEqual(m_particles.size(), 0);

        m_diffusionType = type;

    }

    void swapMainSolverEventWith(KMCEvent *event);

private:

    DiffusionTypes m_diffusionType = DiffusionTypes::TST;

    double m_targetConcentration;

    MainLattice *m_mainLattice;

    Site**** sites;

    uint m_boundaryPadding;

    uint m_NX;
    uint m_NY;
    uint m_NZ;

    uvec3 m_N;


    vector<SoluteParticle*> m_particles;

    double m_kTot;

    vector<double> m_accuAllRates;

    vector<Reaction*> m_allPossibleReactions;

    vector<uint>   m_availableReactionSlots;



    void checkRefCounter();

    void checkAllRefCounters();

    void onConstruct();

    void finalizeObject();


    static bool m_dumpXYZ;
    static bool m_dumpLAMMPS;

    DumpFile *m_dumpFileEvent;
    lammpswriter *m_lammpswriter;

    string m_filepath = "outfiles/";

    SolverEvent *m_solverEvent;

    static uint refCounter;



};


}
