TEMPLATE = app
CONFIG += console

INCLUDEPATH += $$PWD/include

LIBS += -L$$shadowed($$PWD)/lib -lkMC

