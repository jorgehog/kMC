TEMPLATE = subdirs
CONFIG += ordered

SUBDIRS += \
           Quasi2DLoaded \
           WLMC \
#           tests \
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
