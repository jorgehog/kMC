#pragma once

#define _KMCDebugger_IGNORE(expr) static_cast<void>(expr)

//Ignore everything if we are not debugging.
#define KMCDebugger_GetReaction(which) \
    _KMCDebugger_IGNORE(0)
#define KMCDebugger_SetFilename(filename) \
    _KMCDebugger_IGNORE(filename)
#define KMCDebugger_SetPath(path) \
    _KMCDebugger_IGNORE(path)
#define KMCDebugger_SetActiveReaction(reaction) \
    _KMCDebugger_IGNORE(reaction)
#define KMCDebugger_Init() \
    _KMCDebugger_IGNORE(0)
#define KMCDebugger_Finalize() \
    _KMCDebugger_IGNORE(0)
#define KMCDebugger_DumpTrace(i) \
    _KMCDebugger_IGNORE(i)
#define KMCDebugger_DumpFullTrace(...) \
    _KMCDebugger_IGNORE(0)
#define KMCDebugger_MarkPartialStep(_msg) \
    _KMCDebugger_IGNORE(_msg)
#define KMCDebugger_SearchReactionTraceBefore(i) \
    _KMCDebugger_IGNORE(i)
#define KMCDebugger_SearchReactionTraceAfter(i) \
    _KMCDebugger_IGNORE(i)
#define KMCDebugger_SearchImplicationTrace(i) \
    _KMCDebugger_IGNORE(i)
#define KMCDebugger_PushTraces() \
    _KMCDebugger_IGNORE(0)
#define KMCDebugger_PushImplication(site, _pre, _new) \
    _KMCDebugger_IGNORE(site); \
    _KMCDebugger_IGNORE(_pre); \
    _KMCDebugger_IGNORE(_new)
