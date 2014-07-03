include(../app_defaults.pri)

TARGET  = WLMC

SOURCES = WLMCmain.cpp \
    WLMCEvent.cpp \
    wlmcwindow.cpp


OTHER_FILES += infiles/WLMC.cfg

HEADERS += \
    WLMCEvent.h \
    wlmcwindow.h
