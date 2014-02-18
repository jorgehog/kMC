include(../defaults.pri)

TARGET = kMC-tests

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -lunittest++ -L$$TOP_OUT_PWD/src/libs -lkMC

SOURCES = testmain.cpp \
    testbed.cpp \
    snapshot.cpp

HEADERS += \
    testbed.h \
    snapshot.h

OTHER_FILES += \
    ../infiles/knowncase.cfg
