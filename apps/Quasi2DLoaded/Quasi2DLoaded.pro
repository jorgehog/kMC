include(../app_defaults.pri)

TARGET  = Quasi2DLoaded

SOURCES = Quasi2DLoadedmain.cpp \
    ConcentrationControl/concentrationcontrol.cpp

OTHER_FILES += infiles/Quasi2DLoaded.cfg

HEADERS += \
    quasidiffusion.h \
    quasidiffusionevents.h \
    ConcentrationControl/concentrationcontrol.h
