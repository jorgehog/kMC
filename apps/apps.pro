TEMPLATE = subdirs
CONFIG += ordered

SUBDIRS += \
           Quasi2DLoaded \
           tests \
           centerCrystal \
           surfaceGrowth \
           diamondSquareSurface \
           ignisKMC \
           layerGrowth \
           movingSteppes \
           impureSurface \
           stressedSurface \
           WLMC


OTHER_FILES += defaults/default.pro.bones \
               defaults/defaultmain.cpp.bones \
               defaults/defaultconfig.cfg.bones
