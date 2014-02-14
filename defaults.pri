LIBS += -llapack -larmadillo -lconfig++
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += RNG_ZIG

COMMON_CXXFLAGS = -std=c++0x
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_CXXFLAGS -g
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS -O3 -DNDEBUG -DARMA_NO_DEBUG


# Directories
INCLUDEPATH += $$TOP_PWD/src/libs $(HOME)/Dropbox/libs
SRC_DIR = $$TOP_PWD


DIRS = outfiles

for(DIR, DIRS) {
     mkcommands += $$OUT_PWD/$$DIR
}

copydata.commands = $(COPY_DIR) $$PWD/infiles $$OUT_PWD
createDirs.commands = $(MKDIR) $$mkcommands

first.depends = $(first) copydata createDirs
export(first.depends)
export(copydata.commands)
export(createDirs.commands)

QMAKE_EXTRA_TARGETS += first copydata createDirs


CONFIG(RNG_ZIG) {
    DEFINES += KMC_RNG_ZIG
}
