include(../app_defaults.pri)

TARGET  = diamondSquareSurface

SOURCES = diamondSquareSurfacemain.cpp

OTHER_FILES += infiles/diamondSquareSurface.cfg

### DiamondSquare interface

INCLUDEPATH += $$UTILS/diamondSquare

HEADERS += $$UTILS/diamondSquare/src/random/random.h \
           $$UTILS/diamondSquare/src/diamondSquare/diamondSquare.h



SOURCES += $$UTILS/diamondSquare/src/random/random.cpp \
           $$UTILS/diamondSquare/src/diamondSquare/diamondSquare.cpp




