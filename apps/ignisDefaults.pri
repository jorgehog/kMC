INCLUDEPATH += $$PWD/../ignis/include $$PWD/events

LIBS += -L$$PWD/../ignis/lib -lignis

DEFINES += NO_DCVIZ

HEADERS += events/solverevent.h
