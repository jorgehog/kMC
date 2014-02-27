#pragma once

#ifndef KMC_NO_DEBUG

#include "../kmcdebugger.h"

#include "kmcdebugger_class.h"

#include "intrinsicmacros.h"

#include <assert.h>
#include <sstream>

//MISC
#define KMCDebugger_Init(_solverObject) \
    KMCDebugger::initialize(_solverObject)

#define KMCDebugger_Finalize \
    KMCDebugger::reset


#define KMCDebugger_SetFilename(filename) \
    KMCDebugger::traceFileName = filename

#define KMCDebugger_SetPath(path) \
    KMCDebugger::traceFilePath = path


#define KMCDebugger_Assert(A, OP, B, ...) \
    ((A OP B) \
    ? _KMCDebugger_IGNORE(0) \
    : KMCDebugger::_assert(A, B, #OP, #A, #B, \
    __FILE__, \
    __PRETTY_FUNCTION__, \
    __LINE__, ##__VA_ARGS__))

#define KMCDebugger_AssertBool(expr, ...) \
    KMCDebugger_Assert(expr, ==, true, ##__VA_ARGS__)

#define KMCDebugger_AssertClose(A, B, lim, ...) \
    ((A > B) \
    ? KMCDebugger_Assert(A - B, <=, lim, ##__VA_ARGS__) \
    : KMCDebugger_Assert(B - A, <=, lim, ##__VA_ARGS__))


#define KMCDebugger_SetEnabledTo(state) \
    KMCDebugger::setEnabledTo(state)


//TRACE OUTPUT/FETCH FUNCTIONS
#define KMCDebugger_SearchReactionTraceBefore(i) _KMCDebugger_TRACE_SEARCH(KMCDebugger::reactionTraceBefore, i)
#define KMCDebugger_SearchReactionTraceAfter(i)  _KMCDebugger_TRACE_SEARCH(KMCDebugger::reactionTraceAfter, i)
#define KMCDebugger_SearchImplicationTrace(i)    _KMCDebugger_TRACE_SEARCH(KMCDebugger::implicationTrace, i)


//TRACE DUMP CALLERS
#define KMCDebugger_DumpFullTrace(...) \
    KMCDebugger::dumpFullTrace(__LINE__, __FILE__, ##__VA_ARGS__)

#define KMCDebugger_DumpTrace(i) \
    KMCDebugger::dumpPartialTrace(i)


//DEBUGGER FEEDER FUNCTIONS
#define KMCDebugger_SetActiveReaction(reaction) \
    KMCDebugger::setActiveReaction(reaction)

#define KMCDebugger_MarkPartialStep(_msg) \
    KMCDebugger::markPartialStep(_msg)

#define KMCDebugger_PushImplication(site, _pre, _new) \
    KMCDebugger::pushImplication(site, _pre, _new)

#define KMCDebugger_PushTraces() \
    KMCDebugger::pushTraces()

#endif
