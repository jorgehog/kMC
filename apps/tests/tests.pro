include(../app_defaults.pri)

TARGET = tests

LIBS += -lunittest++

SOURCES = testsmain.cpp \
    testbed.cpp \
    snapshot.cpp


HEADERS += testbed.h \
    snapshot.h \
    defines.h

OTHER_FILES += \
    ../../infiles/knowncase.cfg
