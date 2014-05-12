include(../defaults.pri)

TEMPLATE = app
CONFIG += console

INCLUDEPATH  += $$TOP_PWD/include

QMAKE_LIBDIR += $$TOP_PWD/lib

LIBS += -lkMC

DIRS = outfiles

for(DIR, DIRS) {
     mkcommands += $$OUT_PWD/$$DIR
}


QMAKE_EXTRA_TARGETS += first

appDir = $$section(OUT_PWD, "/", 0, -2)

!equals(PWD, $${appDir}) {

    splitdir = $$split(OUT_PWD, "/")
    appname  = $$last(splitdir)

    copydata.commands = $(COPY_DIR) $$PWD/$$appname/infiles $$OUT_PWD

    first.depends += copydata

    QMAKE_EXTRA_TARGETS += copydata
}

createDirs.commands = $(MKDIR) $$mkcommands

first.depends = $(first) $$first.depends createDirs

QMAKE_EXTRA_TARGETS += createDirs
