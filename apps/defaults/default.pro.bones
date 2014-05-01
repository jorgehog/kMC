include(../app_defaults.pri)

TARGET  = __name__

SOURCES = __name__main.cpp


OTHER_FILES += infiles/__name__.cfg

!equals(PWD, $${OUT_PWD}) {

    copydata.commands = $(COPY_DIR) $$PWD/infiles $$OUT_PWD

    QMAKE_EXTRA_TARGETS += copydata
}

createDirs.commands = $(MKDIR) $$mkcommands

QMAKE_EXTRA_TARGETS += createDirs
