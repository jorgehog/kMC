#pragma once

#include <vector>
#include <string>
#include <set>
#include <sys/types.h>

#include <exception>

#include <armadillo>
using arma::wall_clock;
using namespace std;


namespace kMC
{


class Reaction;
class Site;
class SoluteParticle;

class Debugger
{
public:

    static bool enabled;
    static bool prevState;

    static vector<string> reactionTraceBefore;
    static vector<string> reactionTraceAfter;
    static vector<string> implicationTrace;
    static vector<double>      timerData;
    static string implications;

    static uint traceCount;
    static uint implicationCount;

    static Reaction * currentReaction;
    static Reaction * lastCurrentReaction;

    static string reactionString;

    static string _pre;

    static string traceFileName;
    static string traceFilePath;

    static wall_clock timer;

    static set<SoluteParticle*, function<bool(SoluteParticle*, SoluteParticle*)> > affectedUnion;


    //CALLED FROM MACROS
    static void setFilename(const string &filename);
    static void setFilepath(const string &filepath);
    static void setEnabledTo(bool state);
    static void resetEnabled();
    static void pushTraces();
    static void pushImplication(const SoluteParticle *particle, const char *_new);
    static void markPartialStep(const char * msg);
    static void setActiveReaction(Reaction * reaction);
    static void initialize();
    static void reset();
    static void popAffected(SoluteParticle *particle);
    static string fullTrace(int line, const string filename, const string additionalInfo = "");
    static string partialTrace(const uint & i);
    //

    template<typename T>
    static string stringify(const T & val)
    {
        stringstream s;
        s << val;

        return s.str();
    }

    template<typename T>
    static void queuePre(const T & val)
    {
        _pre = stringify(val);
    }

    static string setupAffectedUnion();

    static string addFlagsToImplications();

    static void dumpFullTrace(int line, const char *filename, const string additionalInfo = "");
    static void dumpPartialTrace(const int &i);
};

}
