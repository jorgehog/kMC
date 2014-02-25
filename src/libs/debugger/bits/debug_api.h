#pragma once

#include "../kmcdebugger.h"

#include "kmcdebugger_class.h"

#include "intrinsicmacros.h"

#include <assert.h>
#include <sstream>

//MISC
#define KMCDebugger_Init(_solverObject) \
    KMCDebugger::solverObject = _solverObject; \
    KMCDebugger::timer.tic()

#define KMCDebugger_Finalize KMCDebugger::reset
#define KMCDebugger_SetFilename(filename) KMCDebugger::traceFileName = filename
#define KMCDebugger_SetPath(path)         KMCDebugger::traceFilePath = path
#define KMCDebugger_Assert(A, OP, B, ...) \
    ((A OP B) \
    ? _KMCDebugger_IGNORE(0) \
    : KMCDebugger::_assert(A, B, #OP, #A, #B, \
    __FILE__, \
    __PRETTY_FUNCTION__, \
    __LINE__, ##__VA_ARGS__))
#define KMCDebugger_AssertBool(expr, ...) \
    KMCDebugger_Assert(expr, ==, true, ##__VA_ARGS__)
//

//TRACE OUTPUT/FETCH FUNCTIONS
#define KMCDebugger_SearchReactionTraceBefore(i) _KMCDebugger_TRACE_SEARCH(KMCDebugger::reactionTraceBefore, i)
#define KMCDebugger_SearchReactionTraceAfter(i)  _KMCDebugger_TRACE_SEARCH(KMCDebugger::reactionTraceAfter, i)
#define KMCDebugger_SearchImplicationTrace(i)    _KMCDebugger_TRACE_SEARCH(KMCDebugger::implicationTrace, i)
#define KMCDebugger_GetReaction(which)            KMCDebugger::which##Reaction

#define KMCDebugger_DumpTrace(i) \
    assert(i >= 0); \
    KMCDebugger::dumpPartialTrace(i)

#define KMCDebugger_DumpFullTrace(...) \
    if (KMCDebugger::currentReaction != NULL) {KMCDebugger_PushTraces();} \
    KMCDebugger::dumpFullTrace(__LINE__, __FILE__, ##__VA_ARGS__)
//

//DEBUGGER FEEDER FUNCTIONS
#define KMCDebugger_SetActiveReaction(reaction) \
    KMCDebugger::currentReaction = reaction; \
    KMCDebugger::lastCurrentReaction = reaction; \
    KMCDebugger::reactionString = _KMCDebugger_REACTION_STR(); \
    KMCDebugger::affectedUnion.clear()

#define KMCDebugger_MarkPartialStep(_msg) \
    _KMCDebugger_MAKE_SEPARATOR(_msg); \
    \
    KMCDebugger::implications += KMCDebugger::s.str(); \
    KMCDebugger::implicationCount = 0; \
    \
    _KMCDebugger_CLEAN_SS()

#define KMCDebugger_PushImplication(site, _pre, _new) \
    KMCDebugger::s << "[Implication " << KMCDebugger::implicationCount << "]" << std::endl; \
    KMCDebugger::s << _KMCDebugger_MAKE_IMPLICATION_MESSAGE(site, _pre, _new); \
    \
    KMCDebugger::s << "New affected site(s):\n"; \
    KMCDebugger::setupAffectedUnion(); \
    \
    KMCDebugger::s << "[End of implication " << KMCDebugger::implicationCount << "]\n" << std::endl; \
    KMCDebugger::implications += KMCDebugger::s.str(); \
    \
    KMCDebugger::currentReaction = NULL; \
    KMCDebugger::implicationCount++; \
    \
    _KMCDebugger_CLEAN_SS()

#define KMCDebugger_PushTraces() \
    KMCDebugger::t = KMCDebugger::timer.toc(); \
    KMCDebugger::timerData.push_back(KMCDebugger::t); \
    KMCDebugger::implicationTrace.push_back(KMCDebugger::implications); \
    KMCDebugger::reactionTraceBefore.push_back(KMCDebugger::reactionString); \
    \
    (KMCDebugger::currentReaction != NULL) \
    ? KMCDebugger::reactionTraceAfter.push_back(_KMCDebugger_SITE_STR(KMCDebugger::currentReaction->reactionSite())) \
    : KMCDebugger::reactionTraceAfter.push_back(""); \
    \
    KMCDebugger::implications = _KMCDebugger_INITIAL_IMPLICATION_MSG; \
    KMCDebugger::implicationCount = 0; \
    KMCDebugger::reactionString = "No Reaction Selected"; \
    KMCDebugger::traceCount++; \
    KMCDebugger::timer.tic()
//
