include(../defaults.pri)

TEMPLATE = app
CONFIG += console

INCLUDEPATH  += $$TOP_PWD/include /usr/local/hdf5/include $$PWD

LIBS += -L$$shadowed($$TOP_PWD)/lib -lkMC \
        -L$$UTILS/DCViz/lib -lDCViz \
        -L$$UTILS/HDF5Wrapper/lib -lHDF5Wrapper \
        -L/usr/local/hdf5/lib -lhdf5_cpp -lhdf5

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
