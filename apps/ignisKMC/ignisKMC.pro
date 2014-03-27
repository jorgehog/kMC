include(../app_defaults.pri)
include(../ignisDefaults.pri)

TARGET  = ignisKMC

SOURCES = ignisKMCmain.cpp


OTHER_FILES += infiles/ignisKMC.cfg


copydata.commands = $(COPY_DIR) $$PWD/infiles $$OUT_PWD
createDirs.commands = $(MKDIR) $$mkcommands

first.depends = $(first) copydata createDirs
export(first.depends)
export(copydata.commands)
export(createDirs.commands)

QMAKE_EXTRA_TARGETS += first copydata createDirs
