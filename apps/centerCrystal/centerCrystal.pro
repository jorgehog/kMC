include(../app_defaults.pri)

TARGET  = centerCrystal

SOURCES = centerCrystalmain.cpp


OTHER_FILES += infiles/centerCrystal.cfg


copydata.commands = $(COPY_DIR) $$PWD/infiles $$OUT_PWD
createDirs.commands = $(MKDIR) $$mkcommands

first.depends = $(first) copydata createDirs
export(first.depends)
export(copydata.commands)
export(createDirs.commands)


!equals(PWD, $${OUT_PWD}) {
    QMAKE_EXTRA_TARGETS += first copydata
}

QMAKE_EXTRA_TARGETS += createDirs


