!nompi {

    QMAKE_CXX = mpicxx

    QMAKE_LINK = $$QMAKE_CXX

    QMAKE_LFLAGS += $$system(mpicxx --showme:link)
    COMMON_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK

}

include(../app_defaults.pri)

TARGET  = WLMC

SOURCES = WLMCmain.cpp \
    kmcwlmcsystem.cpp

OTHER_FILES += infiles/WLMC.cfg

HEADERS = kmcwlmcsystem.h

INCLUDEPATH += /usr/local/hdf5/include

LIBS += -L$$UTILS/WLMC/lib -lWLMC \
    -L$$UTILS/HDF5Wrapper/lib -lHDF5Wrapper -lboost_mpi
