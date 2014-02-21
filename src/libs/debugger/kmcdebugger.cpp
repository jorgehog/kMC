#include "kmcdebugger.h"

#ifndef KMC_NO_DEBUG

std::vector<std::string> KMCDebugger::reactionTrace;
std::vector<std::string> KMCDebugger::implicationTrace;
std::string              KMCDebugger::implications = _KMCDebugger_INITIAL_IMPLICATION_MSG;

uint KMCDebugger::traceCount = 0;
uint KMCDebugger::implicationCount = 0;

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

    cout << "---[" << i << " / " << traceCount-1 << " ]" << endl;
    cout << reactionTrace.at(i) << endl;
    cout << implicationTrace.at(i) << endl;

}

#endif
