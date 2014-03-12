#pragma once

#define _KMCDebugger_IGNORE(expr) static_cast<void>(expr)

//Ignore everything if we are not debugging.
#define KMCDebugger_MarkPre(pre) \
    _KMCDebugger_IGNORE(pre)
#define KMCDebugger_SetEnabledTo(state) \
    _KMCDebugger_IGNORE(state)
#define KMCDebugger_IsEnabled \
    false
#define KMCDebugger_Assert(A, OP, B, ...) \
    _KMCDebugger_IGNORE(0)
#define KMCDebugger_AssertBool(expr, ...) \
    _KMCDebugger_IGNORE(expr)
#define KMCDebugger_AssertBreak(...) \
    _KMCDebugger_IGNORE(0)
#define KMCDebugger_AssertClose(A, B, lim, ...) \
    _KMCDebugger_IGNORE(0)
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
#define KMCDebugger_PushImplication(site, _new) \
    _KMCDebugger_IGNORE(site); \
    _KMCDebugger_IGNORE(_new)
