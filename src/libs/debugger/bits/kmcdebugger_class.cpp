#include "kmcdebugger_class.h"

#ifndef KMC_NO_DEBUG

#include <fstream>
#include <sys/time.h>

#include "intrinsicmacros.h"

std::vector<std::string> KMCDebugger::reactionTraceBefore;
std::vector<std::string> KMCDebugger::reactionTraceAfter;
std::vector<std::string> KMCDebugger::implicationTrace;
std::vector<double>      KMCDebugger::timerData;

std::string KMCDebugger::implications = _KMCDebugger_INITIAL_IMPLICATION_MSG;
std::string KMCDebugger::reactionString = _KMCDebugger_INITIAL_REACTION_STR;

Reaction* KMCDebugger::currentReaction = NULL;

uint KMCDebugger::traceCount = 0;
uint KMCDebugger::implicationCount = 0;

std::string KMCDebugger::traceFileName = "";
std::string KMCDebugger::traceFilePath = "";

wall_clock KMCDebugger::timer;

std::stringstream KMCDebugger::s;

double KMCDebugger::t;


void KMCDebugger::dumpFullTrace(bool toFile)
{
    using namespace std;

    if (toFile)
    {
        ofstream file;

        stringstream path;

        if (!traceFilePath.empty())
        {
            path << traceFilePath << "/";
        }

        path << "trace_";

        if (!traceFileName.empty())
        {
            path << traceFileName;
        }

        else
        {
            path << time(NULL);
        }

        path << ".txt";

        file.open(path.str().c_str());
        file << fullTrace() << endl;
        file.close();

        return;

    }

    cout << fullTrace() << endl;

}

void KMCDebugger::dumpPartialTrace(const uint &i)
{
    using namespace std;

    cout << partialTrace(i) << endl;
}

string KMCDebugger::fullTrace()
{
    using namespace std;

    stringstream s;

    s << "======================\n TRACING REACTIONS \n==================\n" << endl;


    for (uint i = 0; i < traceCount; ++i)
    {
        s << partialTrace(i);
    }


    s << "======================\n   TRACE FINISHED   \n==================\n" << endl;

    return s.str();

}

string KMCDebugger::partialTrace(const uint &i)
{
    using namespace std;

    stringstream s;

    s << "---[" << i << " / " << traceCount-1 << " ] " << timerData.at(i)*1000 << " ms" << endl;
    s << reactionTraceBefore.at(i) << endl;
    s << implicationTrace.at(i);

    if (!reactionTraceAfter.at(i).empty())
    {
        s << "\nEnd of reaction look from initial reaction site view: " << endl;
        s << reactionTraceAfter.at(i);
    }
    return s.str();

}

void KMCDebugger::reset()
{

    reactionTraceBefore.clear();
    reactionTraceAfter.clear();
    implicationTrace.clear();
    timerData.clear();

    currentReaction = NULL;

    implications = _KMCDebugger_INITIAL_IMPLICATION_MSG;
    reactionString = _KMCDebugger_INITIAL_REACTION_STR;

    traceCount = 0;
    implicationCount = 0;

    (void)timer.toc();

}

#endif
