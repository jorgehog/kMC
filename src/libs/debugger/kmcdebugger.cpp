#include "kmcdebugger.h"

#include <fstream>
#include <sys/time.h>

#ifndef KMC_NO_DEBUG

std::vector<std::string> KMCDebugger::reactionTrace;
std::vector<std::string> KMCDebugger::implicationTrace;
std::vector<double>      KMCDebugger::timerData;
std::string              KMCDebugger::implications = _KMCDebugger_INITIAL_IMPLICATION_MSG;

uint KMCDebugger::traceCount = 0;
uint KMCDebugger::implicationCount = 0;

wall_clock KMCDebugger::timer;



void KMCDebugger::dumpFullTrace(bool toFile)
{
    using namespace std;

    if (toFile)
    {
        ofstream file;

        stringstream path;

#ifdef KMCDebugger_Path
        path << KMCDebugger_Path << "/";
#endif

        path << "trace" << time(NULL) << ".txt";

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
    s << reactionTrace.at(i) << endl;
    s << implicationTrace.at(i);

    return s.str();

}

void KMCDebugger::reset()
{

    reactionTrace.clear();
    implicationTrace.clear();
    timerData.clear();

    traceCount = 0;
    implicationCount = 0;

    implications = _KMCDebugger_INITIAL_IMPLICATION_MSG;

    (void)timer.toc();

}

#endif
