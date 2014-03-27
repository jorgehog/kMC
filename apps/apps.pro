TEMPLATE = subdirs
CONFIG += ordered

SUBDIRS += tests \
           centerCrystal \
           surfaceGrowth \
           diamondSquareSurface \
           realChalkSetup \
           ignisKMC #__next_app__

OTHER_FILES += defaults/default.pro.bones \
               defaults/defaultmain.cpp.bones \
               defaults/defaultconfig.cfg.bones \
               events/solverevent.h
