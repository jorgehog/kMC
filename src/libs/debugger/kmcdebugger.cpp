#include "kmcdebugger.h"

#ifndef KMC_NO_DEBUG

std::vector<std::string> KMCDebugger::reactionTrace;
std::vector<std::string> KMCDebugger::implicationTrace;
std::vector<double>      KMCDebugger::timerData;
std::string              KMCDebugger::implications;

uint KMCDebugger::traceCount;
uint KMCDebugger::implicationCount;

wall_clock KMCDebugger::timer;

void KMCDebugger::dumpFullTrace()
{
    using namespace std;

    cout << "======================\n TRACING REACTIONS \n==================\n" << endl;


    for (uint i = 0; i < traceCount; ++i)
    {
        dumpPartialTrace(i);
    }


    cout << "======================\n   TRACE FINISHED   \n==================\n" << endl;

}

void KMCDebugger::dumpPartialTrace(const uint &i)
{
    using namespace std;

    cout << "---[" << i << " / " << traceCount-1 << " ] " << timerData.at(i)*1000 << " ms" << endl;
    cout << reactionTrace.at(i) << endl;
    cout << implicationTrace.at(i) << endl;

}

void KMCDebugger::reset()
{
    timer.toc();

    reactionTrace.clear();
    implicationTrace.clear();
    timerData.clear();

    traceCount = 0;
    implicationCount = 0;

    implications = _KMCDebugger_INITIAL_IMPLICATION_MSG;
}

#endif
