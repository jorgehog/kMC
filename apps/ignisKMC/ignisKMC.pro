include(../app_defaults.pri)

TARGET  = ignisKMC

SOURCES = ignisKMCmain.cpp


OTHER_FILES += infiles/ignisKMC.cfg

first.depends = $(first)

!equals(PWD, $${OUT_PWD}) {

    copydata.commands = $(COPY_DIR) $$PWD/infiles $$OUT_PWD

    first.depends += copydata

    QMAKE_EXTRA_TARGETS += copydata
}

createDirs.commands = $(MKDIR) $$mkcommands

first.depends += createDirs

QMAKE_EXTRA_TARGETS += createDirs
INCLUDEPATH += /usr/include/python2.7

LIBS += -lpython2.7
