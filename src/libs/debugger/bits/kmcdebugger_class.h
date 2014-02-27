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

    static KMCSolver * solverObject;

    static std::vector<std::string> reactionTraceBefore;
    static std::vector<std::string> reactionTraceAfter;
    static std::vector<std::string> implicationTrace;
    static std::vector<double>      timerData;
    static std::string implications;

    static uint traceCount;
    static uint implicationCount;

    static std::set<Site*> affectedUnion;

    static void setupAffectedUnion()
    {
        using namespace std;

        vector<Site*> intersect;

        bool currentSiteInPrev;


        for (Site * currentSite : Site::affectedSites())
        {

            currentSiteInPrev = false;
            for (Site * previousSite : affectedUnion)
            {
                if (currentSite == previousSite)
                {
                    currentSiteInPrev = true;
                    break;
                }
            }

            if (!currentSiteInPrev)
            {
                intersect.push_back(currentSite);
            }
        }

        int X, Y, Z;
        for (Site * site : intersect)
        {
            s << "   -" << site->str();

            if (currentReaction != NULL)
            {
                site->distanceTo(currentReaction->reactionSite(), X, Y, Z);
                s << " [" << X << ", " << Y  << ", " << Z << "]";
            }

            s << "\n";

            affectedUnion.insert(site);
        }

        s << "Total: " << intersect.size() << endl;
        intersect.clear();

    }

    static void addFlagsToImplications()
    {
#ifdef KMC_VERBOSE_DEBUG
        stringstream ss, sitess;

        for (Site * site : affectedUnion)
        {
            bool addSite = false;

            sitess << "   " << site->str() << "\n";
            for (Reaction * r : site->siteReactions())
            {
                //we are not interested in logging blocked reactions
                if (!r->isNotBlocked())
                {
                    continue;
                }

                sitess << "    .    " << r->str() << " Flags: ";

                for (int flag : r->updateFlags())
                {
                    if (flag != Reaction::defaultUpdateFlag)
                    {
                        sitess << flag << " ";
                        addSite = true;
                    }
                }

                sitess << "\n";

            }

            if (addSite)
            {
                ss << sitess.str();
                sitess.str(string());
            }
        }

        if (ss.str().empty())
        {
            return;
        }

        implications += ("New non-default flags given to affected site reactions:\n\n" + ss.str());

#endif
    }

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

        cerr << _cerr.str();

        dumpFullTrace(line, file, _cerr.str(), true);

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
