include(../app_defaults.pri)

TARGET  = WLMC

SOURCES = WLMCmain.cpp \
    kmcwlmcsystem.cpp

OTHER_FILES += infiles/WLMC.cfg

HEADERS = kmcwlmcsystem.h

INCLUDEPATH += /usr/local/hdf5/include

LIBS += -L$$UTILS/WLMC/lib -lWLMC \
        -L$$UTILS/HDF5Wrapper/lib -lHDF5Wrapper
LIBS +=-L/usr/local/hdf5/lib -lhdf5_cpp -lhdf5
