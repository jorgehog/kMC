include(../../defaults.pri) 

TEMPLATE = lib 
TARGET = ../../lib/kMC

HEADERS = RNG/kMCRNG.h \
    reactions/reaction.h \
    kmcsolver.h \
    site.h \
    reactions/diffusion/diffusionreaction.h \
    debugger/bits/nodebug.h \
    debugger/bits/intrinsicmacros.h \
    debugger/bits/debug_api.h \
    debugger/debugger.h \
    debugger/bits/debugger_class.h \
    boundary/boundary.h \
    boundary/periodic/periodic.h \
    boundary/concentrationwall/concentrationwall.h \
    boundary/edge/edge.h \
    boundary/wall/surface.h

SOURCES += \
    reactions/reaction.cpp \
    kmcsolver.cpp \
    site.cpp \
    reactions/diffusion/diffusionreaction.cpp \
    RNG/kMCRNG.cpp \
    debugger/bits/debugger_class.cpp \
    boundary/boundary.cpp \
    boundary/periodic/periodic.cpp \
    boundary/concentrationwall/concentrationwall.cpp \
    boundary/wall/surface.cpp

RNG_ZIG {

HEADERS += RNG/zigrandom.h \
           RNG/zignor.h

SOURCES += RNG/zigrandom.cpp \
           RNG/zignor.cpp

}
