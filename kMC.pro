TEMPLATE = subdirs
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += ordered

SUBDIRS += ignis src apps

OTHER_FILES += infiles/config.cfg \
               include/kMC \
               include/kMC_ignis \
               .gitignore
