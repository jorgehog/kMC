include(../app_defaults.pri)

TARGET  = layerGrowth

SOURCES = layerGrowthmain.cpp


OTHER_FILES += infiles/layerGrowth.cfg


copydata.commands = $(COPY_DIR) $$PWD/infiles $$OUT_PWD
createDirs.commands = $(MKDIR) $$mkcommands

first.depends = $(first) copydata createDirs
export(first.depends)
export(copydata.commands)
export(createDirs.commands)

QMAKE_EXTRA_TARGETS += first copydata createDirs
