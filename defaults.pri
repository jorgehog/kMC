LIBS += -llapack -larmadillo -lconfig++
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += RNG_ZIG

QMAKE_LFLAGS_DEBUG += -g

COMMON_CXXFLAGS = -std=c++0x

QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_CXXFLAGS -DKMC_VERBOSE_DEBUG
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS -O3 -DNDEBUG -DARMA_NO_DEBUG -DKMC_NO_DEBUG

INCLUDEPATH += $(HOME)/Dropbox/libs

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
