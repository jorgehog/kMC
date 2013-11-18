include(../../defaults.pri) 

TEMPLATE = lib 
TARGET = kMC

HEADERS = solver.h \
          RNG/kMCRNG.h

SOURCES += solver.cpp

RNG_ZIG {

HEADERS += RNG/zigrandom.h \
           RNG/zignor.h

SOURCES += RNG/zigrandom.cpp \
           RNG/zignor.cpp

}

release {
    QMAKE_CXXFLAGS -= -O2 -O1
    QMAKE_CXXFLAGS += -O3
    DEFINES += ARMA_NO_DEBUG
}
