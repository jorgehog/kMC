include(../../defaults.pri) 

TEMPLATE = lib 
TARGET = kMC

HEADERS = RNG/kMCRNG.h \
    reactions/reaction.h \
    kmcsolver.h \
    site.h \
    reactions/diffusion/diffusionreaction.h

SOURCES += \
    reactions/reaction.cpp \
    kmcsolver.cpp \
    site.cpp \
    reactions/diffusion/diffusionreaction.cpp

RNG_ZIG {

HEADERS += RNG/zigrandom.h \
           RNG/zignor.h

SOURCES += RNG/zigrandom.cpp \
           RNG/zignor.cpp

}


