include(../app_defaults.pri)

TARGET  = movingSteppes

SOURCES = movingSteppesmain.cpp


OTHER_FILES += infiles/movingSteppes.cfg


!equals(PWD, $${OUT_PWD}) {

    copydata.commands = $(COPY_DIR) $$PWD/infiles $$OUT_PWD

    QMAKE_EXTRA_TARGETS += copydata
}

createDirs.commands = $(MKDIR) $$mkcommands

QMAKE_EXTRA_TARGETS += createDirs
