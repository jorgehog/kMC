include(../../defaults.pri) 

TEMPLATE = app 
SOURCES = kmcmain.cpp

LIBS += -L$$TOP_OUT_PWD/src/libs -lkMC
