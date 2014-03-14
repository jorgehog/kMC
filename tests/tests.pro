include(../defaults.pri)
include(../app_defaults.pri)

TARGET = kMC-tests

LIBS += -lunittest++

SOURCES = testmain.cpp \
    testbed.cpp \
    snapshot.cpp

HEADERS += \
    testbed.h \
    snapshot.h \
    defines.h

OTHER_FILES += \
    ../infiles/knowncase.cfg
