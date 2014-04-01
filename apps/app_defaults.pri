include(../defaults.pri)

TEMPLATE = app
CONFIG += console

INCLUDEPATH  += \
    $$TOP_PWD/include \
    $$TOP_PWD/ignis/include

QMAKE_LIBDIR += \
    $$TOP_PWD/lib \
    $$TOP_PWD/ignis/lib

LIBS += \
    -lkMC \
    -lignis


DIRS = outfiles

for(DIR, DIRS) {
     mkcommands += $$OUT_PWD/$$DIR
}
