__author__ = 'jorgehog'

import sys

from matplotlib.pylab import *

sys.path.append("..")
from parse_data import ParseKMCHDF5


obj = ParseKMCHDF5(sys.argv[1])

cases = []

for potential, eb, beta, es, em, h0, c, ignis_index_map, data, n in obj:

    if h0 not in cases:
        cases.append(h0)


for h0_c in cases:

    all_beta_eb = []
    all_c_eq = []
    all_sizes = []

    for potential, eb, beta, es, em, h0, c, ignis_index_map, data, n in obj:

        if h0 != h0_c:
            continue

        c_eq = data.attrs["eqConc"]
        s = data["ignisData"][ignis_index_map["SurfaceSize"], -1]

        if (c_eq < 0):
            print "Warning: negative concentration", c_eq, beta, eb

        all_c_eq.append(c_eq)
        all_sizes.append(s)
        all_beta_eb.append(beta*eb)

    figure(1)
    all_beta_eb, all_sizes, all_c_eq = zip(*[[beb, s, c] for beb, s, c in sorted(zip(all_beta_eb, all_sizes, all_c_eq),
                                                                                 key=lambda x: x[0])])
    plot(all_beta_eb, all_sizes, '-*', label="h0=%g" % h0_c)

    figure(2)
    plot(all_beta_eb, all_c_eq, '-x', label="h0=%g" % h0_c)

figure(1)
xlabel("beta*eb")
ylabel("<s>")
legend()
savefig("sizes_betaeb.png")

save("all_beta_eb", all_beta_eb)
save("all_sizes", all_sizes)
save("all_eqc", all_c_eq)

figure(2)
plot(all_beta_eb, np.exp(-2*array(all_beta_eb)))
xlabel("Eb/kT")
ylabel("C_0/l_0")
legend()
savefig("eqc_betaeb.png")

show()


