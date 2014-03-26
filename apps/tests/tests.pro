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



copydata.commands = $(COPY_DIR) $$PWD/infiles $$OUT_PWD
createDirs.commands = $(MKDIR) $$mkcommands

first.depends = $(first) copydata createDirs
export(first.depends)
export(copydata.commands)
export(createDirs.commands)

QMAKE_EXTRA_TARGETS += first copydata createDirs
