include(../../defaults.pri) 

TEMPLATE = lib 
TARGET = kMC

HEADERS = solver.h \
          RNG/kMCRNG

SOURCES += solver.cpp


RNG_ZIG {

HEADERS += RNG/zigrandom.h \
           RNG/zignor.h

SOURCES += RNG/zigrandom.cpp \
           RNG/zignor.cpp

}
