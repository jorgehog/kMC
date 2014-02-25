#pragma once

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
    static std::vector<std::string> reactionTraceBefore;
    static std::vector<std::string> reactionTraceAfter;
    static std::vector<std::string> implicationTrace;
    static std::vector<double>      timerData;
    static std::string implications;

    static uint traceCount;
    static uint implicationCount;

    static void dumpFullTrace(int line, const char *filename, const string additionalInfo = "", bool toFile = false);
    static void dumpPartialTrace(const uint & i);

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
            std::string what = "")
    {
        using namespace std;

        stringstream s;
        string replString;

        s  << " " << OP << " " << B;

        replString = s.str();
        searchRepl(replString, " == true", "");

        cerr <<  file << ":" << line << ":\n" << func << ":\n";

        cerr << "Assertion '" << A << replString << "' failed: ";

        s.str(string());

        s << "!" << OP;

        replString = s.str();

        searchRepl(replString, "!==", "!=");
        searchRepl(replString, "!!=", "==");
        searchRepl(replString, "!>" , "<=");
        searchRepl(replString, "!<" , ">=");
        searchRepl(replString, "!>=",  "<");
        searchRepl(replString, "!<=",  ">");


        cerr << Aval << " " <<  replString << " " << Bval << ".";

        if (!what.empty())
        {
            cerr << "\nwhat? : " << what;
        }

        cerr << endl;

        exit(1);

    }

    static std::string fullTrace(int line, const string filename, const string additionalInfo = "");
    static std::string partialTrace(const uint & i);

    static void reset();

    static Reaction * currentReaction;
    static Reaction * lastCurrentReaction;

    static std::string reactionString;

    static std::string traceFileName;
    static std::string traceFilePath;

    static wall_clock timer;

    static std::stringstream s;

    static double t;

};
