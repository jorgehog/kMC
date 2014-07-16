include(../app_defaults.pri)

TARGET  = WLMC

SOURCES = WLMCmain.cpp \
    WLMCEvent.cpp \
    kmcwlmcsystem.cpp \
    system.cpp \
    window.cpp


OTHER_FILES += infiles/WLMC.cfg

HEADERS += \
    WLMCEvent.h \
    kmcwlmcsystem.h \
    window.h \
    system.h \
    windowparams.h
