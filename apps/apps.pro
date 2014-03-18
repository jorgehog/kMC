TEMPLATE = subdirs
CONFIG += ordered

SUBDIRS += tests \
           centerCrystal \
           diamondSquareSurface \
           surfaceGrowth \
           realChalkSetup #__next_app__

OTHER_FILES += defaults/default.pro.bones \
               defaults/defaultmain.cpp.bones \
               defaults/defaultconfig.cfg.bones
