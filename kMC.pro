TEMPLATE = subdirs
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += ordered

SUBDIRS += ignis src apps

OTHER_FILES += $$PWD/infiles/config.cfg \
               $$PWD/include/kMC \
               $$PWD/.gitignore
