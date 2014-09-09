import h5py
import sys
import re
from matplotlib.pylab import *

f = h5py.File(sys.argv[1], 'r')

alldata = {}

c0 = 0.9
em0 = 2

n = len(f.keys())
k = 0
for l, run in f.items():
    for potential, data in run.items():

        eb, t, em, c = [float(re.findall("%s\_(\d+\.?\d*)" % ID, potential)[0]) for ID in ["Eb", "beta", "EsMax", "concentration"]]

        if c != c0 or em != em0*eb:
            continue

        if t not in alldata.keys():
            alldata[t] = {}
        if em not in alldata[t].keys():
            alldata[t][em] = {}
        if c not in alldata[t][em].keys():
            alldata[t][em][c] = {}

        alldata[t][em][c][l] = {}

        for i, name in enumerate(data["ignisEventDescriptions"][0]):
            alldata[t][em][c][l][str(name).split("@")[0]] = data["ignisData"][i, :]
        alldata[t][em][c][l]["name"] = potential

    print k+1, " / ", n, "complete."
    k += 1


for t, rest in alldata.items():
    for em, rest2 in rest.items():
        for c, lengths in rest2.items():
                figure()

                hold("on")
                for length, data in lengths.items():
                    loglog(data["Time"], data["heightRMS"], label=length)

                    title(data["name"])

                legend(loc=2)
                xlabel("t [s * R_0]")
                ylabel("RMS(h)")
                draw()
                savefig(data["name"] + ".png")

f.close()

