include(../app_defaults.pri)

TARGET  = WLMC

SOURCES = WLMCmain.cpp \
    WLMCEvent.cpp \
    wlmcwindow.cpp \
    wlmcsystem.cpp


OTHER_FILES += infiles/WLMC.cfg

HEADERS += \
    WLMCEvent.h \
    wlmcwindow.h \
    wlmcsystem.h
