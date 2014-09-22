import h5py
from hgext.mq import clone
import sys
import re
from gst._gst import TIME_ARGS
from matplotlib.pylab import *
from math import log
from numpy import where, asarray, random



def interpolate(value_array, time_array, time, start=0):

    if start >= len(time_array):
        return start, 0, 0

    i = start
    while time_array[i] < time:
        i += 1

        if i == len(time_array):
            return i, 0, 0

    if i == 0:

        if time != time_array[0]:
            return 0, 0, 0

        return i, 1, value_array[0]

    incline = (value_array[i] - value_array[i-1])/(time_array[i] - time_array[i-1])

    return i, 1, value_array[i-1] + incline*(time - time_array[i-1])


def combine_results(time_array, value_array, do_roughen=True):

    if do_roughen:
        time_array = roughen(*time_array)
        value_array = roughen(*value_array)

    times = set()
    for time in time_array:
        time = asarray(time)
        times = times.union(set(time))

    print min(times), max(times)

    times = array(sorted(list(times)))
    times = roughen(times)

    all_f_t = zeros([len(value_array), len(times)])
    scales = zeros(len(times))
    starts = zeros(len(value_array))

    length = 40
    for ti, time in enumerate(times):
        for i in range(len(value_array)):



            loc, scale, all_f_t[i][ti] = interpolate(value_array[i], time_array[i], time, starts[i])

            scales[ti] += scale
            starts[i] = loc

        if scales[ti] == 0:
            print "ERROR IN SCALING", ti, time
            sys.exit(1)

        s = (ti*length)/(len(times)-1)
        print "\rcombining: [%s%s]" % ("-"*s, " "*(length-s)),
        sys.stdout.flush()

    return times, all_f_t.sum(0)/scales


def roughen(*args):

    n_args = len(args)

    new = []
    for i in range(n_args):

        idx = [0]
        k = 1
        N = len(args[i])

        while k < N:
            idx.append(int(round(k)))
            k += 0.001*k*log(1 + k)

        idx = unique(asarray(idx, dtype=int))

        new.append(args[i][idx])

    if n_args == 1:
        return new[0]
    else:
        return tuple(new)



def get_w_beta(rms, times, do_plot=False):

    t, f_t = combine_results(times, rms, do_roughen=True)

    from scipy import polyfit

    M = 10000000000000000
    M2 = M

    logt = np.log(t)
    logf = np.log(f_t)
    N = len(logt)

    I = None
    I2 = None
    R = None
    R2 = None
    for i in range(N/10, N):
        r, residuals, c, c, c = polyfit(logt[:i], logf[:i], 1, full=True)

        if residuals[0] < M:
            M = residuals[0]
            R = r
            I = i

        i2 = N - i
        r2, residuals2, c, c, c = polyfit(logt[i2:-1], logf[i2:-1], 1, full=True)

        if residuals2[0] < M2:
            M2 = residuals2[0]
            R2 = r2
            I2 = i2

    if do_plot:
        cla()
        loglog(t, f_t)

        figure()
        plot(logt, logf)
        plot(logt[:I], R[1] + logt[:I]*R[0])
        plot(logt[I2:], R2[1] + logt[I2:]*R2[0])

        print "w inf = %g, beta = %g" % (R[0], R2[1])

        show()

    return R[0], R2[1]


def main():

    f = h5py.File(sys.argv[1], 'r')

    all_data = {}

    c0 = 0.4
    em0 = 10
    t0 = 0.1
    n0 = 1
    l0 = 1024

    total = len(f.keys())
    k = 0
    for l, run in f.items():

        print "l=%5s %2d / %2d Parsing %d items ..." % (l, k + 1, total, len(run.items()))

        for potential, data in run.items():

            if not "concentration" in potential:
                potential += "_concentration_0.4"

            eb, t, em, c = [float(re.findall("%s\_(\d+\.?\d*)" % ID, potential)[0]) for ID in ["Eb",
                                                                                               "beta",
                                                                                               "EsMax",
                                                                                               "concentration"]]

            n = re.findall("\_n_(\d+)", potential)

            if n:
                n = int(n[0])
            else:
                n = ""

            # print c,c0, em, em0*eb, t, t0, l, l0
            if c != c0 or em != em0*eb or t != t0:
                # if c != c0:
                #     print "c"
                # if em != em0*eb:
                #     print "em"
                # if t != t0:
                #     print "t"
                # if l != l0:
                #     print "l"
                continue

            if t not in all_data.keys():
                all_data[t] = {}
            if em not in all_data[t].keys():
                all_data[t][em] = {}
            if c not in all_data[t][em].keys():
                all_data[t][em][c] = {}
            if l not in all_data[t][em][c].keys():
                all_data[t][em][c][l] = {}

            all_data[t][em][c][l][n] = {}

            for i, name in enumerate(data["ignisEventDescriptions"][0]):
                all_data[t][em][c][l][n][str(name).split("@")[0]] = data["ignisData"][i, :]
            all_data[t][em][c][l][n]["name"] = potential

        k += 1

    f.close()

    all_times = []
    all_lengths = []
    all_rms = []

    for t, rest in all_data.items():
        for em, rest2 in rest.items():
            for c, rest3 in rest2.items():

                figure()
                hold("on")

                name = rest3.values()[0].values()[0]["name"].split("_n_")[0]

                print name,
                for length, all_runs in rest3.items():
                    print length

                    time_array = [data["Time"] for data in all_runs.values()]
                    value_array = [data["heightRMS"] for data in all_runs.values()]

                    toss = []
                    for i, t in enumerate(time_array):
                        if where(t < 0)[0].size != 0:
                            print "Negative times in run ", i
                            toss.append(i)

                    for i in toss[::-1]:
                        time_array.pop(i)
                        value_array.pop(i)




                    t, f_t = combine_results(time_array, value_array, do_roughen=True)

                    loglog(t, f_t, label=length)
                    all_times.append(t)
                    all_lengths.append(length)
                    all_rms.append(f_t)

                print
                title(name)

                legend(loc=2)
                xlabel("t [s * R_0]")
                ylabel("<RMS(h)>")
                draw()
                savefig(name + ".png")

                clf()

                # winf, beta = get_w_beta(all_rms, all_times)

    #
    # close("all")
    # for t, rest in all_data.items():
    #     for em, rest2 in rest.items():
    #         for c, lengths in rest2.items():
    #             figure()
    #
    #             name = lengths.values()[0]["name"]
    #
    #             hold("on")
    #             for length, data in lengths.items():
    #                 plot(data["Time"], data["height"]/data["Time"], label=length)
    #
    #             xscale("log")
    #             title(name)
    #
    #             legend(loc=2)
    #             xlabel("t [s * R_0]")
    #             ylabel("<v(t)>")
    #             draw()
    #
    #
    # show()




if __name__ == "__main__":
    main()