CONFIG -= app_bundle
CONFIG -= qt
CONFIG += RNG_ZIG


QMAKE_CXX = gcc

!noccache {
    QMAKE_CXX = ccache $$QMAKE_CXX
}

COMMON_CXXFLAGS = -std=c++11 -fopenmp

QMAKE_CXXFLAGS += \
    $$COMMON_CXXFLAGS

QMAKE_CXXFLAGS_DEBUG += \
    $$COMMON_CXXFLAGS \
    -DKMC_VERBOSE_DEBUG

QMAKE_CXXFLAGS_RELEASE += \
    $$COMMON_CXXFLAGS \
    -O3 \
    -DNDEBUG \
    -DARMA_NO_DEBUG \
    -DKMC_NO_DEBUG

QMAKE_LFLAGS_RELEASE -= -O1
QMAKE_LFLAGS_RELEASE += -O3

INCLUDEPATH += \
    $(HOME)/Dropbox/libs \
    $$PWD/ignis/include

QMAKE_LIBDIR += $$PWD/ignis/lib

LIBS += \
    -larmadillo \
    -lconfig++ \
    -lignis \
    -fopenmp

DEFINES += \
    ARMA_MAT_PREALLOC=3

CONFIG(NO_OMP) {
    DEFINES += KMC_NO_OMP
}

CONFIG(RNG_ZIG) {

    DEFINES += KMC_RNG_ZIG

}


CONFIG(compphys) {

    INCLUDEPATH  += $(HOME)/shared/code
    INCLUDEPATH  += $(HOME)/shared/armadillo-4.200.0/usr/include
    QMAKE_LIBDIR += $(HOME)/shared/armadillo-4.200.0/usr/lib

}

TOP_PWD = $$PWD
