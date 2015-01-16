import sys
import os

import matplotlib.pylab as plt

from os.path import join

lib_path = join(os.getcwd(), "pyblockzor", "lib")
sys.path.append(lib_path)

import pyblockzor

sys.path.append(join(os.getcwd(), ".."))
from parse_data import ParseKMCHDF5

path = sys.argv[1]

all_results = {}

parser = ParseKMCHDF5(path)

for potential, alpha, mu, E0, s0, r0, ignis_index_map, data, n in parser:

    if not "SurfaceSizeLocal" in ignis_index_map.keys():
        continue

    surface_size = data["ignisData"][ignis_index_map["SurfaceSizeLocal"]]

    print "Blocking %s of size %d" % (potential, len(surface_size))

    satisfied = False
    while not satisfied:

        consistent = False
        while not consistent:

            rawparams = raw_input("minbs maxbs nbs = ")

            minbs, maxbs, nbs = [int(x) for x in rawparams.split()]

            consistent = pyblockzor.consistencyCheck(minbs, maxbs, nbs, len(surface_size), False)

        sizes, sigmas = pyblockzor.block(surface_size, minbs, maxbs, nbs)

        plt.plot(sizes, sigmas)
        plt.title(potential)
        plt.xlabel("nblocks")
        plt.ylabel("sigma")

        plt.show()

        ans = raw_input("satisfied? [y]")

        satisfied = ans == "y"



