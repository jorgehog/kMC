TEMPLATE = subdirs
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += ordered

SUBDIRS += utils src apps

src.depends = utils
apps.depends = src utils

OTHER_FILES += infiles/config.cfg \
               include/kMC \
               .gitignore
