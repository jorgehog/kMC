#pragma once

#include "debugger_class.h"

#include "../../site.h"

#include "../../reactions/reaction.h"

#define _KMCDebugger_INITIAL_IMPLICATION_MSG "[Implications]: \n"
#define _KMCDebugger_INITIAL_REACTION_STR    "Initialization"

#define _KMCDebugger_TRACE_SEARCH(trace, i) \
    ((i < 0) \
    ? trace.at((kMC::Debugger::traceCount + (i))) \
    : trace.at(i))

#ifdef KMC_VERBOSE_DEBUG

#define _KMCDebugger_REACTION_STR() \
    kMC::Debugger::currentReaction->info()

#define _KMCDebugger_SITE_STR(site) \
    site->info()

#else

#define _KMCDebugger_REACTION_STR() \
    kMC::Debugger::currentReaction->str()

#define _KMCDebugger_SITE_STR(site) \
    site->str()

#endif

#define _KMCDebugger_MAKE_IMPLICATION_MESSAGE(site, _pre, _new) \
    (_KMCDebugger_SITE_STR(site) + \
    ("  What? " + \
    ((std::string)_pre + (" -> " + (std::string)_new + \
    ("\n\n")))))
