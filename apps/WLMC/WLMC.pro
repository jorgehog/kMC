include(../app_defaults.pri)

TARGET  = WLMC

SOURCES = WLMCmain.cpp \
    kmcwlmcsystem.cpp

OTHER_FILES += infiles/WLMC.cfg

HEADERS = kmcwlmcsystem.h

LIBS += -L$$UTILS/WLMC/lib -lWLMC

