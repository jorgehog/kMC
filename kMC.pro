TEMPLATE = subdirs
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += ordered

SUBDIRS += src tests

OTHER_FILES += $$PWD/infiles/config.cfg \
               $$PWD/include/kMC \
               $$PWD/.qmake.conf \
               $$PWD/.gitignore

