#pragma once

#include "debug_api.h"

#include <vector>
#include <string>
#include <set>
#include <sys/types.h>

#include <exception>

#include <armadillo>
using arma::wall_clock;

class Reaction;
class Site;

using namespace std;

class KMCDebugger
{
public:

    static bool enabled;

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

    static string traceFileName;
    static string traceFilePath;

    static wall_clock timer;

    static set<Site*> affectedUnion;

    //CALLED FROM MACROS
    static void setFilename(const string &filename);
    static void setFilepath(const string &filepath);
    static void setEnabledTo(bool state);
    static void pushTraces();
    static void pushImplication(Site *site, const char *_pre, const char *_new);
    static void markPartialStep(const char * msg);
    static void setActiveReaction(Reaction * reaction);
    static void initialize();
    static void reset();
    static string fullTrace(int line, const string filename, const string additionalInfo = "");
    static string partialTrace(const uint & i);
    //

    static string setupAffectedUnion();

    static string addFlagsToImplications();

    static void dumpFullTrace(int line, const char *filename, const string additionalInfo = "");
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
            string what = "",
            string additionalInfo = "")
    {

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
        searchRepl(replString, "!>=",  "<");
        searchRepl(replString, "!<=",  ">");
        searchRepl(replString, "!>" , "<=");
        searchRepl(replString, "!<" , ">=");


        _cerr << Aval << " " <<  replString << " " << Bval << ".";

        if (!what.empty())
        {
            _cerr << "\nwhat? : " << what;
        }


        cerr << _cerr.str() << endl;


        dumpFullTrace(line, file, additionalInfo);

    }




};
