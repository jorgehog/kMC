__author__ = 'jorgehog'

import sys

from matplotlib.pylab import *

sys.path.append("..")
from parse_data import ParseKMCHDF5



obj = ParseKMCHDF5(sys.argv[1])

all_E0 = []
for potential, alpha, mu, E0, s0, r0, ignis_index_map, data, n in obj:
    if E0 not in all_E0:
        all_E0.append(E0)

for E0_case in all_E0:

    all_alpha = []
    all_mu_eq = []
    all_sizes = []

    for potential, alpha, mu, E0, s0, r0, ignis_index_map, data, n in obj:

        if E0 != E0_case:
            continue

        if "muEq" in data.attrs:
            mu_eq = data.attrs["muEq"]
        else:
            mu_eq = mu

        s = data["ignisData"][ignis_index_map["SurfaceSize"], -1]

        all_mu_eq.append(mu_eq)
        all_sizes.append(s)
        all_alpha.append(alpha)

    figure(1)
    all_alpha, all_sizes, all_mu_eq = zip(*[[beb, s, c] for beb, s, c in sorted(zip(all_alpha, all_sizes, all_mu_eq),
                                                                                 key=lambda x: x[0])])
    plot(all_alpha, all_sizes, '-*', label="E0=%g" % E0)

    figure(2)
    plot(all_alpha, all_mu_eq, '-x', label="E0=%g" % E0)

figure(1)
xlabel("alpha")
ylabel("<s>")
legend()
savefig("sizes_betaeb.png")

save("all_beta_eb", all_alpha)
save("all_sizes", all_sizes)
save("all_eqc", all_mu_eq)

figure(2)
xlabel("alpha")
ylabel("mu")
legend()
savefig("eqc_betaeb.png")

show()


