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

        cerr <<  file << ":" << line << ":\n" << func << ":\n";

        cerr << "Assertion '" << A << " " << OP << " " << B << "' failed: ";

        cerr << Aval << " !" << OP << " " << Bval << ".";

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
