#pragma once

#ifndef KMC_NO_DEBUG

#include "../debugger.h"

#include "debugger_class.h"

#include "intrinsicmacros.h"

#include <assert.h>
#include <sstream>

//MISC
#define KMCDebugger_Init() \
    kMC::Debugger::initialize()

#define KMCDebugger_Finalize() \
    kMC::Debugger::reset()

#define KMCDebugger_IsEnabled \
    kMC::Debugger::enabled


#define KMCDebugger_SetFilename(filename) \
    kMC::Debugger::setFilename(filename)

#define KMCDebugger_SetPath(path) \
    kMC::Debugger::setFilepath(path)


#define KMCDebugger_Assert(A, OP, B, ...) \
    (((A) OP (B)) \
    ? static_cast<void>(0) \
    : kMC::Debugger::_assert(A, B, #OP, #A, #B, \
    __FILE__, \
    __PRETTY_FUNCTION__, \
    __LINE__, ##__VA_ARGS__))

#define KMCDebugger_AssertClose(A, B, lim, ...) \
    ((A > B) \
    ? KMCDebugger_Assert(A - B, <=, lim, ##__VA_ARGS__) \
    : KMCDebugger_Assert(B - A, <=, lim, ##__VA_ARGS__))

#define KMCDebugger_AssertBool(expr, ...) \
    KMCDebugger_Assert(expr, ==, true, ##__VA_ARGS__)

#define KMCDebugger_AssertBreak(...) \
    KMCDebugger_AssertBool(false, ##__VA_ARGS__)

#define KMCDebugger_SetEnabledTo(state) \
    kMC::Debugger::setEnabledTo(state)


//TRACE OUTPUT/FETCH FUNCTIONS
#define KMCDebugger_SearchReactionTraceBefore(i) _KMCDebugger_TRACE_SEARCH(kMC::Debugger::reactionTraceBefore, i)
#define KMCDebugger_SearchReactionTraceAfter(i)  _KMCDebugger_TRACE_SEARCH(kMC::Debugger::reactionTraceAfter, i)
#define KMCDebugger_SearchImplicationTrace(i)    _KMCDebugger_TRACE_SEARCH(kMC::Debugger::implicationTrace, i)


//TRACE DUMP CALLERS
#define KMCDebugger_DumpFullTrace(...) \
    kMC::Debugger::dumpFullTrace(__LINE__, __FILE__, #__VA_ARGS__)

#define KMCDebugger_DumpTrace(i) \
    kMC::Debugger::dumpPartialTrace(i)


//DEBUGGER FEEDER FUNCTIONS
#define KMCDebugger_MarkPre(pre) \
    kMC::Debugger::queuePre(pre);

#define KMCDebugger_SetActiveReaction(reaction) \
    kMC::Debugger::setActiveReaction(reaction)

#define KMCDebugger_MarkPartialStep(_msg) \
    kMC::Debugger::markPartialStep(_msg)

#define KMCDebugger_PushImplication(site, _new) \
    kMC::Debugger::pushImplication(site,  _new)

#define KMCDebugger_PushTraces() \
    kMC::Debugger::pushTraces()

#endif
