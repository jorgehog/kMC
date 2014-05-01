include(../app_defaults.pri)

TARGET  = ignisKMC

SOURCES = ignisKMCmain.cpp


OTHER_FILES += infiles/ignisKMC.cfg

!equals(PWD, $${OUT_PWD}) {

    copydata.commands = $(COPY_DIR) $$PWD/infiles $$OUT_PWD

    QMAKE_EXTRA_TARGETS += copydata
}

createDirs.commands = $(MKDIR) $$mkcommands

QMAKE_EXTRA_TARGETS += createDirs

INCLUDEPATH += $(HOME)/code/DCViz/include

INCLUDEPATH += /usr/include/python2.7

LIBS += -lpython2.7
