include(../defaults.pri)

TEMPLATE = lib

TARGET = ../lib/kMC


QMAKE_LFLAGS += -g

QMAKE_CXXFLAGS_RELEASE += -g

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
    boundary/periodicshifted/periodicshifted.h

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
    boundary/periodicshifted/periodicshifted.cpp


RNG_ZIG {

HEADERS += RNG/zigrandom.h \
           RNG/zignor.h

SOURCES += RNG/zigrandom.cpp \
           RNG/zignor.cpp

}

!equals(PWD, $${OUT_PWD}) {
    QMAKE_POST_LINK += $(COPY_DIR) $$OUT_PWD/../lib $$TOP_PWD
}



