__author__ = 'jorgehog'

import sys

from matplotlib.pylab import *

sys.path.append("..")
from parse_data import ParseKMCHDF5


obj = ParseKMCHDF5(sys.argv[1])

all_beta_eb = []
all_c_eq = []
all_sizes = []
for potential, eb, beta, em, c, ignis_index_map, data, n in obj:
    c_eq = data.attrs["eqConc"]
    s = data["ignisData"][ignis_index_map["SurfaceSize"], -1]

    all_c_eq.append(c_eq)
    all_sizes.append(s)
    all_beta_eb.append(beta*eb)

figure()
plot(all_beta_eb, all_sizes, 'k*')
xlabel("beta*eb")
ylabel("<s>")
savefig("sizes_betaeb.png")

save("all_beta_eb", all_beta_eb)
save("all_sizes", all_sizes)
save("all_eqc", all_c_eq)

figure()
plot(all_beta_eb, all_c_eq, 'kx')
xlabel("Eb/kT")
ylabel("C_0/l_0")
savefig("eqc_betaeb.png")
show()


