include(../app_defaults.pri)

TARGET  = WLMC

SOURCES = WLMCmain.cpp \
    kmcwlmcsystem.cpp

OTHER_FILES += infiles/WLMC.cfg

HEADERS = kmcwlmcsystem.h

INCLUDEPATH += $$PWD/WLMC/include

LIBS += -L$$PWD/WLMC/lib -lWLMC
