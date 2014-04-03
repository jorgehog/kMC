CONFIG -= app_bundle
CONFIG -= qt
CONFIG += RNG_ZIG

QMAKE_CXX = gcc

COMMON_CXXFLAGS = -std=c++11

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

INCLUDEPATH += \
    $(HOME)/Dropbox/libs \
    $$PWD/ignis/include

QMAKE_LIBDIR += $$PWD/ignis/lib

LIBS += \
    -larmadillo \
    -lconfig++
    -lignis

DEFINES += \
    ARMA_MAT_PREALLOC=3

CONFIG(RNG_ZIG) {

    DEFINES += \
         KMC_RNG_ZIG
}

TOP_PWD = $$PWD
