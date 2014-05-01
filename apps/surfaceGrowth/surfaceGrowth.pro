include(../app_defaults.pri)

TARGET  = surfaceGrowth

SOURCES = surfaceGrowthmain.cpp


OTHER_FILES += infiles/surfaceGrowth.cfg


copydata.commands = $(COPY_DIR) $$PWD/infiles $$OUT_PWD
createDirs.commands = $(MKDIR) $$mkcommands

first.depends = $(first) copydata createDirs
export(first.depends)
export(copydata.commands)
export(createDirs.commands)

QMAKE_EXTRA_TARGETS += first

!equals(PWD, $${OUT_PWD}) {
    QMAKE_EXTRA_TARGETS += copydata
}

QMAKE_EXTRA_TARGETS += createDirs
