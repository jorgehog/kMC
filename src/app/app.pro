include(../../defaults.pri) 

TEMPLATE = app 
SOURCES = kmcmain.cpp

LIBS += -L$$TOP_OUT_PWD/src/libs -lkMC




DIRS = outfiles

for(DIR, DIRS) {
     mkcommands += $$OUT_PWD/$$DIR
}

createDirs.commands = $(MKDIR) $$mkcommands

first.depends += createDirs
QMAKE_EXTRA_TARGETS += createDirs
