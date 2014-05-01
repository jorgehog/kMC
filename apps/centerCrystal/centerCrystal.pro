include(../app_defaults.pri)

TARGET  = centerCrystal

SOURCES = centerCrystalmain.cpp


OTHER_FILES += infiles/centerCrystal.cfg


!equals(PWD, $${OUT_PWD}) {

    copydata.commands = $(COPY_DIR) $$PWD/infiles $$OUT_PWD

    QMAKE_EXTRA_TARGETS += copydata
}

createDirs.commands = $(MKDIR) $$mkcommands

QMAKE_EXTRA_TARGETS += createDirs
