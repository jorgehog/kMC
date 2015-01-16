__author__ = 'jorgehog'

import sys

from scipy.stats import linregress

from matplotlib.pylab import *

sys.path.append("..")
from parse_data import ParseKMCHDF5



obj = ParseKMCHDF5(sys.argv[1])

all_E0 = []
for potential, alpha, mu, E0, s0, r0, ignis_index_map, data, n in obj:
    if E0 not in all_E0:
        all_E0.append(E0)
print all_E0

all_slopes = []

X=0
v=0
shift=[2.5,2, 1.1,1.3,1.4,1.5,1.6]
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
    plot(all_alpha, all_sizes, '-*', label="E0=%g" % E0_case)

    figure(2)
    print "SIGN CORRECTION"
    UBERSAT = exp(all_mu_eq) - 1
    plot(all_alpha, UBERSAT, '-x', label="E0=%g" % E0_case)

    slope, a, b, c, d = linregress(all_alpha, UBERSAT)

    figure(3)
    scatter(abs(E0_case), slope)

    all_slopes.append(slope)
    #
    # figure(4)
    # plot([a*abs(E0_case)**0.15 for a in all_alpha], all_sizes, label="%d E0=%g" % (X, E0_case))
    # X+=1


import numpy
stuff = numpy.loadtxt("/tmp/boltzmann_ascii.arma")


figure(1)
plot(stuff[:, 0], stuff[:, 1], "k--^", label="analytical")
xlabel("alpha")
ylabel("<s>")
legend()
savefig("sizes_betaeb.png")

save("all_beta_eb", all_alpha)
save("all_sizes", all_sizes)
save("all_eqc", all_mu_eq)

figure(2)
xlabel("alpha")
ylabel("supersaturation")
legend(loc=0)
savefig("eqc_betaeb.png")

Slope, a, b, c, d = linregress([-E0 for E0 in all_E0], all_slopes)

figure(3)
xlabel("$|E_0|$")
ylabel(r"$d\sigma/d\alpha$")
title("slope = %g" % Slope)
savefig("supersatslopes.png")

# figure(4)
# legend()

show()


