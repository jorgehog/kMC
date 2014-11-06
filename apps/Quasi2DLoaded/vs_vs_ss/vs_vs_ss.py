__author__ = 'jorgehog'

import sys
import os

from matplotlib.pylab import *

sys.path.append(os.getcwd() + "/..")
from parse_data import ParseKMCHDF5


def plot_c(local_c, c_eq, local_v, beta, eb):
    plot(array(local_c)/c_eq, local_v, "*", markersize=10, label="b Eb=%g %g" % (beta, eb))

obj = ParseKMCHDF5(sys.argv[1])

local_c = []
local_v = []
k = 0

eb_prev = 0
beta_prev = 0
c_eq = 1

figure()
for potential, eb, beta, em, c, ignis_index_map, data, n in obj:
    if data.attrs["useConcEquil"] == 1:
        c_eq = data.attrs["eqConc"]
        print c_eq, local_c
        continue

    if (beta != beta_prev or eb != eb_prev) and len(local_c) != 0:

        print k, beta, beta_prev, eb, eb_prev, len(local_c) != 0
        plot_c(local_c, c_eq, local_v, beta, eb)

        local_c = []
        local_v = []

        k += 1

    beta_prev = beta
    eb_prev = eb

    local_c.append(c)

    h_final = data["ignisData"][ignis_index_map["height"], -1]
    h_start = data["ignisData"][ignis_index_map["height"], 0]
    T = data["ignisData"][ignis_index_map["Time"], -1]

    v = (h_final - h_start)/T

    local_v.append(v)

save("local_c", array(local_c)/c_eq)
save("local_v", local_v)
plot_c(local_c, c_eq, local_v, beta, eb)
xlabel("C/C_0")
ylabel("v_s")
legend()
grid()
show()