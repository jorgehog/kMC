include(../app_defaults.pri)

TARGET  = ignisKMC

SOURCES = ignisKMCmain.cpp


OTHER_FILES += infiles/ignisKMC.cfg

INCLUDEPATH += /usr/include/python2.7

LIBS += -lpython2.7
