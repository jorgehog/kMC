TEMPLATE = subdirs
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += ordered

SUBDIRS += src tests

BUILD_DIR = $$OUT_PWD

OTHER_FILES += $$PWD/infiles/config.cfg
