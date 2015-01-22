from sympy.mpmath.functions.zetazeros import sure_number_block
import sys
import os
import h5py

import numpy as np
import matplotlib.pylab as plt

from os.path import join

lib_path = join(os.getcwd(), "pyblockzor", "lib")
sys.path.append(lib_path)

import pyblockzor

sys.path.append(join(os.getcwd(), ".."))
from parse_data import ParseKMCHDF5

path = sys.argv[1]

def get_names(filepath):
    path, name = os.path.split(filepath)

    filename = join(path, ".".join(name.split(".")[:-1]) + "_blocking.dat")

    return path, name, filename

def store_data(potential, s, res, filepath):

    path, name, filename = get_names(filepath)

    with open(filename, 'a') as f:
        f.write("%s %g %g\n" %(potential, s, res))

def is_in(name, filepath):

    path, name_, filename = get_names(filepath)

    if not os.path.exists(filename):
        return False

    with open(filename, 'r') as f:
        for line in f:
            if len(line.split()) != 3:
                continue

            if line.split()[0] == name:
                return True

    return False

all_results = {}

parser = ParseKMCHDF5(path)

for l, potential_raw, alpha, mu, E0, s0, r0, ignis_index_map, data, n in parser:

    potential = "L_%s_%s" % (l, potential_raw)

    if not "SurfaceSizeLocal" in ignis_index_map.keys():
        continue


    print "Blocking %s" % (potential)

    if is_in(potential, path):
        print "found in set."

    run = raw_input("Run?[y]")


    if run != "y" and run != "":
        continue


    surface_size = data["ignisData"][ignis_index_map["SurfaceSizeLocal"]]
    surface_size = surface_size[np.where(surface_size != 0)]

    plt.plot(surface_size[::100])
    plt.plot(np.cumsum(surface_size[::100])/(np.arange(len(surface_size[::100])) + 1), 'k--')
    plt.show()

    therm = raw_input("thermalize=")
    therm = int(therm)*100
    surface_size = surface_size[therm:]

    print "size: ", len(surface_size)

    final_s_surf = surface_size.mean()
    final_s = data["ignisData"][ignis_index_map["SurfaceSize"], -1]
    print "result: ", final_s, " ", final_s_surf
    print "original stddev: ", np.std(surface_size)/len(surface_size)**0.5

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

        satisfied = ans == "y" or ans == ""

    result = raw_input("Result: ")

    all_results[potential] = (final_s_surf, float(result))

    store_data(potential, final_s_surf, float(result), path)

parser.close()