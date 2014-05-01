include(../app_defaults.pri)

TARGET  = surfaceGrowth

SOURCES = surfaceGrowthmain.cpp


OTHER_FILES += infiles/surfaceGrowth.cfg

!equals(PWD, $${OUT_PWD}) {

    copydata.commands = $(COPY_DIR) $$PWD/infiles $$OUT_PWD

    QMAKE_EXTRA_TARGETS += copydata
}

createDirs.commands = $(MKDIR) $$mkcommands

QMAKE_EXTRA_TARGETS += createDirs

