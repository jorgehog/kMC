#pragma once

#include "debug_api.h"

#include "../../reactions/reaction.h"

#include <vector>
#include <string>
#include <sys/types.h>

#include <exception>

#include <armadillo>

using arma::wall_clock;

class KMCDebugger
{
public:

    static bool enabled;

    static KMCSolver * solverObject;

    static std::vector<std::string> reactionTraceBefore;
    static std::vector<std::string> reactionTraceAfter;
    static std::vector<std::string> implicationTrace;
    static std::vector<double>      timerData;
    static std::string implications;

    static uint traceCount;
    static uint implicationCount;

    static Reaction * currentReaction;
    static Reaction * lastCurrentReaction;

    static std::string reactionString;

    static std::string traceFileName;
    static std::string traceFilePath;

    static wall_clock timer;

    static std::set<Site*> affectedUnion;

    //CALLED FROM MACROS
    static void setEnabledTo(bool state);
    static void pushTraces();
    static void pushImplication(Site *site, const char *_pre, const char *_new);
    static void markPartialStep(const char * msg);
    static void setActiveReaction(Reaction * reaction);
    static void initialize(KMCSolver * solver);
    static void reset();
    static std::string fullTrace(int line, const string filename, const string what, const string additionalInfo = "");
    static std::string partialTrace(const uint & i);
    //

    static string setupAffectedUnion();

    static string addFlagsToImplications();

    static void dumpFullTrace(int line, const char *filename, const string what = "", const string additionalInfo = "", bool toFile = true);
    static void dumpPartialTrace(const int &i);

    static void searchRepl(string & s, string _find, string _repl)
    {

        int position = s.find(_find);
        while (position != (int)string::npos)
        {
            s.replace(position, _find.size(), _repl);
            position = s.find(_find, position + 1);
        }

    }


    template<typename TA, typename TB>
    static void
    _assert(TA Aval,
            TB Bval,
            const char * OP,
            const char * A,
            const char * B,
            const char * file,
            const char * func,
            int line,
            std::string what = "",
            std::string additionalInfo = "")
    {
        using namespace std;

        stringstream s, _cerr;
        string replString;

        s  << " " << OP << " " << B;

        replString = s.str();
        searchRepl(replString, " == true", "");

        _cerr <<  file << ":" << line << ":\n" << func << ":\n";

        _cerr << "Assertion '" << A << replString << "' failed: ";

        s.str(string());

        s << "!" << OP;

        replString = s.str();

        searchRepl(replString, "!==", "!=");
        searchRepl(replString, "!!=", "==");
        searchRepl(replString, "!>" , "<=");
        searchRepl(replString, "!<" , ">=");
        searchRepl(replString, "!>=",  "<");
        searchRepl(replString, "!<=",  ">");


        _cerr << Aval << " " <<  replString << " " << Bval << ".";

        if (!what.empty())
        {
            _cerr << "\nwhat? : " << what;
        }

        _cerr << endl;

        if (enabled)
        {
            cerr << _cerr.str();
        }

        dumpFullTrace(line, file, _cerr.str(), additionalInfo);

    }




};
