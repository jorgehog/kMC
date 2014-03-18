include(../defaults.pri)

TEMPLATE = app
CONFIG += console

INCLUDEPATH += $$TOP_PWD/include

LIBS += -L$$TOP_PWD/lib -lkMC

