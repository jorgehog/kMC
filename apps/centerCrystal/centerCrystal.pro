include(../app_defaults.pri)

TARGET  = centerCrystal

SOURCES = centerCrystalmain.cpp


OTHER_FILES += infiles/centerCrystal.cfg


first.depends = $(first)

!equals(PWD, $${OUT_PWD}) {

    copydata.commands = $(COPY_DIR) $$PWD/infiles $$OUT_PWD

    first.depends += copydata

    QMAKE_EXTRA_TARGETS += copydata
}

createDirs.commands = $(MKDIR) $$mkcommands

first.depends += createDirs

QMAKE_EXTRA_TARGETS += createDirs
