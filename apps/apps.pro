TEMPLATE = subdirs
CONFIG += ordered

SUBDIRS += tests \
           centerCrystal \
           surfaceGrowth \
           diamondSquareSurface \
           ignisKMC \
           layerGrowth #__next_app__

OTHER_FILES += defaults/default.pro.bones \
               defaults/defaultmain.cpp.bones \
               defaults/defaultconfig.cfg.bones \
               events/solverevent.h
