TEMPLATE = subdirs
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += ordered

SUBDIRS += ignis \
    WLMC \
    diamondSquare \
    HDF5Wrapper \
    DCViz

HEADERS += lammpswriter/lammpswriter.h \
    libconfig_utils/libconfig_utils.h \
    BADAss/badass.h
