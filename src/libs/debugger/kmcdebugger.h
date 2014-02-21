#pragma once

#ifndef KMC_NO_DEBUG

#include "../reactions/reaction.h"

#include <vector>
#include <string>
#include <sstream>
#include <sys/types.h>

#include <assert.h>

#include <armadillo>

using arma::wall_clock;

class KMCDebugger
{
public:
    static std::vector<std::string> reactionTrace;
    static std::vector<std::string> implicationTrace;
    static std::vector<double>      timerData;
    static std::string implications;

    static uint traceCount;
    static uint implicationCount;

    static void dumpFullTrace();
    static void dumpPartialTrace(const uint & i);
    static void reset();

    static wall_clock timer;
};

#define _KMCDebugger_INITIAL_IMPLICATION_MSG "[implications]: \n"

#define KMCDebugger_Init() \
    KMCDebugger::reset(); \
    KMCDebugger::timer.tic()

#define KMCDebugger_dumpTrace(i) \
    assert(i >= 0); \
    KMCDebugger::dumpPartialTrace(i)

#define KMCDebugger_dumpFullTrace() \
    KMCDebugger::dumpFullTrace()

#define _KMCDebugger_TRACE_SEARCH(trace, i) \
    ((i < 0) \
    ? expr.at((KMCDebugger::traceCount + (i))) \
    : expr.at(i))

#define KMCDebugger_PushInitialTrace() \
    _KMCDebugger_PUSHTRACES_FROMSTR("Initialization")

#ifdef KMC_VERBOSE_DEBUG
#define KMCDebugger_pushTraces(reaction) \
    _KMCDebugger_PUSHTRACES_FROMSTR(reaction->info())
#else
#define KMCDebugger_pushTraces(reaction) \
    _KMCDebugger_PUSHTRACES_FROMSTR(reaction->str())
#endif

#define _KMCDebugger_PUSHTRACES_FROMSTR(_reactionStr) \
    KMCDebugger::timerData.push_back(KMCDebugger::timer.toc()); \
    KMCDebugger::traceCount++; \
    KMCDebugger::implicationTrace.push_back(KMCDebugger::implications); \
    KMCDebugger::implicationCount = 0; \
    KMCDebugger::reactionTrace.push_back(_reactionStr); \
    KMCDebugger::implications = _KMCDebugger_INITIAL_IMPLICATION_MSG; \
    KMCDebugger::timer.tic()

#ifdef KMC_VERBOSE_DEBUG
#define _KMCDebugger_MAKE_IMPLICATION_MESSAGE(site, _pre, _new) \
    (site->info() + ("\nParticleStateChange: " + ((std::string)_pre + ("->" + (std::string)_new + (";\n")))))
#else
#define _KMCDebugger_MAKE_IMPLICATION_MESSAGE(site, _pre, _new) \
    ("  "  + (site->str() + (": " + ((std::string)_pre + ("->" + (std::string)_new + (";\n"))))))
#endif

#define KMCDebugger_pushImplication(site, _pre, _new) \
    KMCDebugger::implications += _KMCDebugger_MAKE_IMPLICATION_MESSAGE(site, _pre, _new); \
    KMCDebugger::implicationCount++

#define _KMCDebugger_MAKE_SEPARATOR(s, _msg) \
    s << "### " << _msg << "  " << "prev. imp.: " << KMCDebugger::implicationCount << " ###\n" << std::endl


#define KMCDebugger_addImplicationSeparator(_msg) \
    std::stringstream s; \
    _KMCDebugger_MAKE_SEPARATOR(s, _msg); \
    KMCDebugger::implications += s.str(); \
    KMCDebugger::implicationCount = 0


#else

#define _KMCDebugger_IGNORE(expr) static_cast<void>(expr)

//Ignore everything if we are not debugging.
#define KMCDebugger_Init() \
    _KMCDebugger_IGNORE(0)
#define KMCDebugger_dumpTrace(i) \
    _KMCDebugger_IGNORE(i)
#define KMCDebugger_dumpFullTrace() \
    _KMCDebugger_IGNORE(0)
#define KMCDebugger_addImplicationSeparator(_msg) \
    _KMCDebugger_IGNORE(_msg)
#define _KMCDebugger_TRACE_SEARCH(trace, i) \
    _KMCDebugger_IGNORE(trace); \
    _KMCDebugger_IGNORE(i)
#define KMCDebugger_pushTraces(reaction) \
    _KMCDebugger_IGNORE(reaction)
#define KMCDebugger_pushImplication(site, _pre, _new) \
    _KMCDebugger_IGNORE(site); \
    _KMCDebugger_IGNORE(_pre); \
    _KMCDebugger_IGNORE(_new)
#define KMCDebugger_PushInitialTrace() \
    _KMCDebugger_IGNORE(0)
#endif


#define KMCDebugger_searchReactionTrace(i) _KMCDebugger_TRACE_SEARCH(KMCDebugger::reactionTrace, i)
#define KMCDebugger_searchImplicationTrace(i) _KMCDebugger_TRACE_SEARCH(KMCDebugger::implicationTrace, i)
