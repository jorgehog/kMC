include(../app_defaults.pri)

TARGET = tests

LIBS += -lunittest++

SOURCES = testsmain.cpp \
    testbed/testbed.cpp \
    snapshot/snapshot.cpp


HEADERS += testbed/testbed.h \
    snapshot/snapshot.h \
    defines.h \
    dummyreaction.h

OTHER_FILES += infiles/knowncase.cfg
