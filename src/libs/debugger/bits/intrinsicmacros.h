#pragma once

#include "kmcdebugger_class.h"
#include "../../reactions/reaction.h"
#include "../../reactions/diffusion/diffusionreaction.h"

#define _KMCDebugger_INITIAL_IMPLICATION_MSG "[implications]: \n"
#define _KMCDebugger_INITIAL_REACTION_STR    "Initialization"

#define _KMCDebugger_TRACE_SEARCH(trace, i) \
    ((i < 0) \
    ? trace.at((KMCDebugger::traceCount + (i))) \
    : trace.at(i))

#define _KMCDebugger_REACTIONCAST() \
    ((KMCDebugger::currentReaction->name == "DiffusionReaction") \
    ? ((DiffusionReaction*)KMCDebugger::currentReaction) \
    : ((Reaction*)KMCDebugger::currentReaction))

#ifdef KMC_VERBOSE_DEBUG
#define _KMCDebugger_REACTION_STR() _KMCDebugger_REACTIONCAST()->info()
#define _KMCDebugger_SITE_STR(site) site->info()
#else
#define _KMCDebugger_REACTION_STR() _KMCDebugger_REACTIONCAST()->str()
#define _KMCDebugger_SITE_STR(site) site->str()
#endif

#define _KMCDebugger_MAKE_IMPLICATION_MESSAGE(site, _pre, _new) \
("  "  + (_KMCDebugger_SITE_STR(site) + (": " + ((std::string)_pre + ("->" + (std::string)_new + (";\n"))))))

#define _KMCDebugger_MAKE_SEPARATOR(_msg) \
    KMCDebugger::s << "##### " << _msg << "  " << "prev. imp.: " << KMCDebugger::implicationCount << " #####\n\n" << std::endl

#define _KMCDebugger_CLEAN_SS() KMCDebugger::s.str(std::string());
