include(../app_defaults.pri)

TARGET  = __name__

SOURCES = __name__main.cpp


OTHER_FILES += infiles/__name__.cfg


copydata.commands = $(COPY_DIR) $$PWD/infiles $$OUT_PWD
createDirs.commands = $(MKDIR) $$mkcommands

first.depends = $(first) copydata createDirs
export(first.depends)
export(copydata.commands)
export(createDirs.commands)


!equals(PWD, $${OUT_PWD}) {
    QMAKE_EXTRA_TARGETS += first copydata createDirs
}
