TEMPLATE = subdirs
CONFIG += ordered

SUBDIRS += tests \
           centerCrystal \
           surfaceGrowth \
           diamondSquareSurface \
           ignisKMC \
           layerGrowth \
           movingSteppes


OTHER_FILES += defaults/default.pro.bones \
               defaults/defaultmain.cpp.bones \
               defaults/defaultconfig.cfg.bones
