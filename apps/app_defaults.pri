include(../defaults.pri)

TEMPLATE = app
CONFIG += console

INCLUDEPATH  += $$TOP_PWD/include

QMAKE_LIBDIR += $$TOP_PWD/lib

LIBS += -lkMC
