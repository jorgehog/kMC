TEMPLATE = subdirs
CONFIG += ordered

SUBDIRS += tests \
           centerCrystal \
           surfaceGrowth \
           diamondSquareSurface \
           ignisKMC \
           layerGrowth \
           movingSteppes \
           impureSurface \
           stressedSurface \
           WLMC \
           Quasi2DLoaded


OTHER_FILES += defaults/default.pro.bones \
               defaults/defaultmain.cpp.bones \
               defaults/defaultconfig.cfg.bones
