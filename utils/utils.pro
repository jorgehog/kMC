TEMPLATE = subdirs
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += ordered

SUBDIRS += ignis \
    WLMC \
    diamondSquare

HEADERS += lammpswriter/lammpswriter.h \
    libconfig_utils/libconfig_utils.h \
    BADAss/BADAss.h
