__author__ = 'jorgehog'

import os
import re
import h5py


class ParseKMCHDF5:

    def __init__(self, which_file):
        self.path, self.filename = os.path.split(which_file)

        self.file = h5py.File(which_file, 'r')

    def __iter__(self):

        for l, run in self.file.items():

            for potential, data in run.items():

                alpha, mu, E0, s0, r0 = [float(re.findall("%s\_(-?\d+\.?\d*)" % ID, potential)[0]) for ID in
                                           ["alpha", "mu", "E0", "s0", "r0"]]

                n = re.findall("\_n_(\d+)", potential)

                _ignis_index_map = {}

                for i, name in enumerate(data["ignisEventDescriptions"][0]):
                    name_strip = str(name).split("@")[0]
                    _ignis_index_map[name_strip] = i

                yield l, potential, alpha, mu, E0, s0, r0, _ignis_index_map, data, n

    def getfile(self):
        return self.file

    def close(self):
        self.file.close()

if __name__ == "__main__":
    obj = ParseKMCHDF5("/home/jorgehog/code/build-kMC-Desktop_Qt_5_3_GCC_64bit-Release/apps/Quasi2DLoaded/Quasi2D.h5")

    for potential, alpha, mu, em, c, ignis_index_map, data, n in obj:
        print alpha, data["ignisData"][ignis_index_map["SurfaceSize"], -1]