include(../app_defaults.pri)

TARGET  = diamondSquareSurface

SOURCES = diamondSquareSurfacemain.cpp

OTHER_FILES += infiles/diamondSquareSurface.cfg


copydata.commands = $(COPY_DIR) $$PWD/infiles $$OUT_PWD
createDirs.commands = $(MKDIR) $$mkcommands

first.depends = $(first) copydata createDirs
export(first.depends)
export(copydata.commands)
export(createDirs.commands)

!equals(PWD, $${OUT_PWD}) {
    QMAKE_EXTRA_TARGETS += first copydata createDirs
}

### DiamondSquare interface

INCLUDEPATH += $$PWD/diamondSquare

HEADERS += diamondSquare/src/random/random.h \
           diamondSquare/src/diamondSquare/diamondSquare.h



SOURCES += diamondSquare/src/random/random.cpp \
           diamondSquare/src/diamondSquare/diamondSquare.cpp




