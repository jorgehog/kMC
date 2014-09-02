TEMPLATE = subdirs
CONFIG += ordered

SUBDIRS += \
           Quasi2DLoaded \
           tests \
#           WLMC \
#           centerCrystal \
#           surfaceGrowth \
#           diamondSquareSurface \
#           ignisKMC \
#           layerGrowth \
#           movingSteppes \
#           impureSurface \
#           stressedSurface \


OTHER_FILES += defaults/default.pro.bones \
               defaults/defaultmain.cpp.bones \
               defaults/defaultconfig.cfg.bones \
               commonkmcevents.h
