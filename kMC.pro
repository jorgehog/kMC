TEMPLATE = subdirs
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += ordered

SUBDIRS += utils src apps

OTHER_FILES += infiles/config.cfg \
               include/kMC \
               .gitignore
