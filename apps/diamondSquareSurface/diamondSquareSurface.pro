include(../app_defaults.pri)

TARGET  = diamondSquareSurface

SOURCES = diamondSquareSurfacemain.cpp

OTHER_FILES += infiles/diamondSquareSurface.cfg


first.depends = $(first)

!equals(PWD, $${OUT_PWD}) {

    copydata.commands = $(COPY_DIR) $$PWD/infiles $$OUT_PWD

    first.depends += copydata

    QMAKE_EXTRA_TARGETS += copydata
}

createDirs.commands = $(MKDIR) $$mkcommands

first.depends += createDirs

QMAKE_EXTRA_TARGETS += createDirs

### DiamondSquare interface

INCLUDEPATH += $$PWD/diamondSquare

HEADERS += diamondSquare/src/random/random.h \
           diamondSquare/src/diamondSquare/diamondSquare.h



SOURCES += diamondSquare/src/random/random.cpp \
           diamondSquare/src/diamondSquare/diamondSquare.cpp




