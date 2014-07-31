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

#define KMCDebugger_SetEnabledTo(state) \
    kMC::Debugger::setEnabledTo(state)

#define KMCDebugger_ResetEnabled() \
    kMC::Debugger::resetEnabled()

//TRACE OUTPUT/FETCH FUNCTIONS
#define KMCDebugger_SearchReactionTraceBefore(i) _KMCDebugger_TRACE_SEARCH(kMC::Debugger::reactionTraceBefore, i)
#define KMCDebugger_SearchReactionTraceAfter(i)  _KMCDebugger_TRACE_SEARCH(kMC::Debugger::reactionTraceAfter, i)
#define KMCDebugger_SearchImplicationTrace(i)    _KMCDebugger_TRACE_SEARCH(kMC::Debugger::implicationTrace, i)

//BADAss interface
#define KMCBAI(info, ...) [&] (const badass::BADAssException &exc) { \
    kMC::Debugger::dumpFullTrace(exc.whichLine(), exc.whichFile().c_str(), kMC::Debugger::stringify(info)); \
    ##__VA_ARGS__ \
    }


//TRACE DUMP CALLERS
#define KMCDebugger_DumpFullTrace(...) \
    kMC::Debugger::dumpFullTrace(__LINE__, __FILE__, #__VA_ARGS__)

#define KMCDebugger_DumpTrace(i) \
    kMC::Debugger::dumpPartialTrace(i)


//DEBUGGER FEEDER FUNCTIONS
#define KMCDebugger_PopAffected(particle) \
    kMC::Debugger::popAffected(particle)

#define KMCDebugger_MarkPre(pre) \
    kMC::Debugger::queuePre(pre)

#define KMCDebugger_SetActiveReaction(reaction) \
    kMC::Debugger::setActiveReaction(reaction)

#define KMCDebugger_MarkPartialStep(_msg) \
    kMC::Debugger::markPartialStep(_msg)

#define KMCDebugger_PushImplication(particle, _new) \
    kMC::Debugger::pushImplication(particle,  _new)

#define KMCDebugger_PushTraces() \
    kMC::Debugger::pushTraces()

#endif
