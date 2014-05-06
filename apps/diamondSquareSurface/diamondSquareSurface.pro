include(../app_defaults.pri)

TARGET  = diamondSquareSurface

SOURCES = diamondSquareSurfacemain.cpp

OTHER_FILES += infiles/diamondSquareSurface.cfg

### DiamondSquare interface

INCLUDEPATH += $$PWD/diamondSquare

HEADERS += diamondSquare/src/random/random.h \
           diamondSquare/src/diamondSquare/diamondSquare.h



SOURCES += diamondSquare/src/random/random.cpp \
           diamondSquare/src/diamondSquare/diamondSquare.cpp




