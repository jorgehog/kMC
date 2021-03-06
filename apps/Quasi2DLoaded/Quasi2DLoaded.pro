include(../app_defaults.pri)

TARGET  = Quasi2DLoaded

SOURCES = Quasi2DLoadedmain.cpp \
    movingwall.cpp \
    quasidiffusionevents.cpp \
    quasidiffusionsystem.cpp

OTHER_FILES += infiles/Quasi2DLoaded.cfg

HEADERS += \
    quasidiffusion.h \
    quasidiffusionevents.h \
    movingwall.h \
    quasidiffusionsystem.h
