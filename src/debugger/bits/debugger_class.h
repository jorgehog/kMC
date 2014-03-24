#pragma once

#include "debug_api.h"

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

    static set<Site*> affectedUnion;

    //CALLED FROM MACROS
    static void setFilename(const string &filename);
    static void setFilepath(const string &filepath);
    static void setEnabledTo(bool state);
    static void resetEnabled();
    static void pushTraces();
    static void pushImplication(Site *site, const char *_new);
    static void markPartialStep(const char * msg);
    static void setActiveReaction(Reaction * reaction);
    static void initialize();
    static void reset();
    static string fullTrace(int line, const string filename, const string additionalInfo = "");
    static string partialTrace(const uint & i);
    //

    template<typename T>
    static void queuePre(const T & val)
    {
        stringstream s;
        s << val;

        _pre = s.str();

    }

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

        stringstream OP_B, assertMessage, assertCore;
        string OP_B_repl, A_repl;


        OP_B  << " " << OP << " " << B;

        OP_B_repl = OP_B.str();
        searchRepl(OP_B_repl, " == true", "");


        A_repl    = string(A);
        searchRepl(A_repl   , "false", "");


        assertMessage <<  file << ":" << line << ": " << func << ": ";

        assertCore << A_repl << OP_B_repl;

        if (!assertCore.str().empty())
        {

            assertMessage << "Assertion '" << assertCore.str() << "' failed: ";

            OP_B.str(string());

            OP_B << "!" << OP;

            OP_B_repl = OP_B.str();

            searchRepl(OP_B_repl, "!==", "!=");
            searchRepl(OP_B_repl, "!!=", "==");
            searchRepl(OP_B_repl, "!>=",  "<");
            searchRepl(OP_B_repl, "!<=",  ">");
            searchRepl(OP_B_repl, "!>" , "<=");
            searchRepl(OP_B_repl, "!<" , ">=");


            assertMessage << Aval << " " <<  OP_B_repl << " " << Bval;

        }

        else
        {
            assertMessage << "Assertion failed";
        }

        if (what.empty())
        {
            assertMessage << ".";
        }
        else
        {
            assertMessage << " : " << what;
        }


        dumpFullTrace(line, file, additionalInfo);

        throw std::runtime_error(assertMessage.str());

    }




};

}
