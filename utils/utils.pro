TEMPLATE = subdirs
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += ordered

SUBDIRS += ignis \
    WLMC \
    diamondSquare \
    HDF5Wrapper

HEADERS += lammpswriter/lammpswriter.h \
    libconfig_utils/libconfig_utils.h \
    DCViz/include/DCViz.h \
    BADAss/BADAss.h
