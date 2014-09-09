import h5py
import sys
import re
from matplotlib.pylab import *

f = h5py.File(sys.argv[1], 'r')

alldata = {}

n = len(f.keys())
k = 0
for l, run in f.items():
    for potential, data in run.items():
        t = float(re.findall("beta\_(\d+\.?\d*)", potential)[0])

        if t not in alldata.keys():
            alldata[t] = {}

        alldata[t][l] = {}

        for i, name in enumerate(data["ignisEventDescriptions"][0]):
            alldata[t][l][str(name).split("@")[0]] = data["ignisData"][i, :]

    print k+1, " / ", n, "complete."
    k += 1

for t, lengths in alldata.items():
    figure()

    title("beta = %g" % t)

    hold("on")
    for length, data in lengths.items():
        loglog(data["Time"], data["heightRMS"], label=length)

    legend(loc=2)
    xlabel("t [s * R_0]")
    ylabel("RMS(h)")
    draw()
    savefig("t_vs_hrms_beta" + str(t) + ".png")

show()

f.close()

