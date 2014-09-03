include(../defaults.pri)

TEMPLATE = lib

TARGET = ../lib/kMC

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
    particlestates.h \
    ignisinterface/solverevent.h \
    ignisinterface/kmcevent.h \
    ignisinterface/kmcparticles.h \
    soluteparticle.h \
    boundary/periodicshifted/periodicshifted.h \
    potential/potential.h \
    potential/stressedsurface/stressedsurface.h \
    potential/inertwall/inertwall.h \
    system.h \
    reactions/diffusion/arrheniusdiffusion.h \
    reactions/diffusion/tstdiffusion.h \
    boundary/edge/sphericalEdge.h

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
    particlestates.cpp \
    soluteparticle.cpp \
    ignisinterface/kmcparticles.cpp \
    ignisinterface/kmcevent.cpp \
    boundary/periodicshifted/periodicshifted.cpp \
    potential/stressedsurface/stressedsurface.cpp \
    potential/inertwall/inertwall.cpp \
    system.cpp \
    reactions/diffusion/arrheniusdiffusion.cpp \
    reactions/diffusion/tstdiffusion.cpp \
    boundary/edge/sphericalEdge.cpp


RNG_ZIG {

HEADERS += RNG/zigrandom.h \
           RNG/zignor.h

SOURCES += RNG/zigrandom.cpp \
           RNG/zignor.cpp

}

QMAKE_PRE_LINK += $(MKDIR) $$PWD/../lib $$shadowed($$PWD)/../lib

!equals(PWD, $${OUT_PWD}) {
    QMAKE_POST_LINK += $(COPY_DIR) $$OUT_PWD/../lib $$TOP_PWD
}



