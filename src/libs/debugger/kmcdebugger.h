#pragma once

#define _KMCDebugger_IGNORE(expr) static_cast<void>(expr)

#ifndef KMC_NO_DEBUG
#include "bits/debug_api.h"
#else
#include "bits/nodebug.h"
#endif

