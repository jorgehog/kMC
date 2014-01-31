#ifndef KMCRNG_H
#define KMCRNG_H

#ifdef KMC_RNG_ZIG

#include "zignor.h"
#include "zigrandom.h"

#define KMC_RNG_NORMAL DRanNormalZig32

#define KMC_RNG_UNIFORM DRan_MWC8222

#define KMC_INIT_RNG(seed)                       \
    int inseed = seed;                       \
    int cseed = 100;                         \
    int seed2 = inseed * 3;                  \
    RanSetSeed_MWC8222(&seed2, cseed);       \
    RanNormalSetSeedZig32(&inseed, 5)       \

#endif

struct Seed {
    enum SeedType {
        fromTime,
        specific
    };
};

#endif //KMCRNG_H
